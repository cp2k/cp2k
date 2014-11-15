!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
MODULE topology_gromos
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                             parser_get_object,&
                                             parser_search_string
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE f77_blas
  USE input_cp2k_restarts_util,        ONLY: section_velocity_val_set
  USE input_section_types,             ONLY: section_get_rval,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE string_table,                    ONLY: s2s,&
                                             str2id
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: atom_info_type,&
                                             connectivity_info_type,&
                                             topology_parameters_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_gromos'

  PRIVATE
  PUBLIC :: read_topology_gromos, read_coordinate_g96

CONTAINS

! *****************************************************************************
!> \brief Read GROMOS topology file
!> \param file_name ...
!> \param topology ...
!> \param para_env ...
!> \param subsys_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE read_topology_gromos (file_name,topology,para_env,subsys_section, error)
    CHARACTER(LEN=*), INTENT(IN)             :: file_name
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_topology_gromos', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: ctemp
    CHARACTER(LEN=default_string_length)     :: label, string
    CHARACTER(LEN=default_string_length), &
      DIMENSION(21)                          :: avail_section
    CHARACTER(LEN=default_string_length), &
      POINTER                                :: namearray1(:), namearray2(:)
    INTEGER :: begin, handle, i, iatom, ibond, icon, ii(50), index_now, &
      iresid, isolvent, itemp, itype, iw, jatom, natom, natom_prev, nbond, &
      nbond_prev, ncon, nsolute, nsolvent, ntype, offset, offset2, stat
    INTEGER, POINTER                         :: ba(:), bb(:), na(:)
    LOGICAL                                  :: failure, found
    REAL(KIND=dp)                            :: ftemp
    REAL(KIND=dp), POINTER                   :: ac(:), am(:)
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_parser_type), POINTER            :: parser

    NULLIFY(parser, logger)
    NULLIFY(namearray1,namearray2)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/GTOP_INFO",&
         extension=".subsysLog",error=error)
    CALL timeset(routineN,handle)

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
    avail_section(21) = "SOLVENTCONSTR"

    atom_info => topology%atom_info
    conn_info => topology%conn_info

    natom_prev = 0
    IF(ASSOCIATED(atom_info%id_molname)) natom_prev = SIZE(atom_info%id_molname)
    ! TITLE SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the TITLE section'
    CALL parser_create(parser,file_name,para_env=para_env,error=error)
    label = TRIM(avail_section(1))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       DO
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,string,string_length=80,error=error)
          IF(string == TRIM("END")) EXIT
          IF(iw>0) WRITE(iw,*) "GTOP_INFO| ",TRIM(string)
       END DO
    END IF
    CALL parser_release(parser,error=error)
    ! TOPPHYSCON SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the TOPPHYSCON section'
    CALL parser_create(parser,file_name,para_env=para_env,error=error)
    label = TRIM(avail_section(2))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       DO
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,string,error=error)
          IF(string == TRIM("END")) EXIT
          IF(iw>0) WRITE(iw,*) "GTOP_INFO|  ",TRIM(string)
       END DO
    END IF
    CALL parser_release(parser,error=error)
    ! TOPVERSION SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the TOPVERSION section'
    CALL parser_create(parser,file_name,para_env=para_env,error=error)
    label = TRIM(avail_section(3))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       DO
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,string,error=error)
          IF(string == TRIM("END")) EXIT
          IF(iw>0) WRITE(iw,*) "GTOP_INFO|  ",TRIM(string)
       END DO
    END IF
    CALL parser_release(parser,error=error)
    ! ATOMTYPENAME SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the ATOMTYPENAME section'
    CALL parser_create(parser,file_name,para_env=para_env,error=error)
    label = TRIM(avail_section(4))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ntype,error=error)
       CALL reallocate(namearray1,1,ntype)
       DO itype=1,ntype
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,namearray1(itype),error=error)
          IF(iw>0) WRITE(iw,*) "GTOP_INFO|  ",TRIM(namearray1(itype))
       END DO
    END IF
    CALL parser_release(parser,error=error)
    ! RESNAME SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the RESNAME section'
    CALL parser_create(parser,file_name,para_env=para_env,error=error)
    label = TRIM(avail_section(5))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ntype,error=error)
       CALL reallocate(namearray2,1,ntype)
       DO itype=1,ntype
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,namearray2(itype),error=error)
          IF(iw>0) WRITE(iw,*) "GTOP_INFO|  ",TRIM(namearray2(itype))
       END DO
    END IF
    CALL parser_release(parser,error=error)
    ! SOLUTEATOM SECTION
    iresid=1
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the SOLUTEATOM section'
    CALL parser_create(parser,file_name,para_env=para_env,error=error)
    label = TRIM(avail_section(6))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,natom,error=error)
       CALL reallocate(atom_info%id_molname,1,natom_prev+natom)
       CALL reallocate(atom_info%resid,1,natom_prev+natom)
       CALL reallocate(atom_info%id_resname,1,natom_prev+natom)
       CALL reallocate(atom_info%id_atmname,1,natom_prev+natom)
       CALL reallocate(atom_info%id_element,1,natom_prev+natom)
       CALL reallocate(atom_info%atm_charge,1,natom_prev+natom)
       CALL reallocate(atom_info%atm_mass,1,natom_prev+natom)
       CALL parser_get_next_line(parser,1,error=error)
       DO iatom=1,natom
          index_now = iatom + natom_prev
          CALL parser_get_object(parser,itemp,error=error)
          CALL parser_get_object(parser,itemp,error=error)
          atom_info%resid(index_now)=itemp
          atom_info%id_molname(index_now)=str2id(s2s(namearray2(itemp)))
          atom_info%id_resname(index_now)=str2id(s2s(namearray2(itemp)))
          CALL parser_get_object(parser,string,error=error)
          CALL parser_get_object(parser,itemp,error=error)
          atom_info%id_atmname(index_now)=str2id(s2s(namearray1(itemp)))
          atom_info%id_element(index_now)=str2id(s2s(namearray1(itemp)))
          CALL parser_get_object(parser,atom_info%atm_mass(index_now),error=error)
          CALL parser_get_object(parser,atom_info%atm_charge(index_now),error=error)
          CALL parser_get_object(parser,itemp,error=error)
          IF(iw>0) WRITE(iw,*) "GTOP_INFO| PUT SOLUTEATOM INFO HERE!!!!"
          CALL parser_get_object(parser,ntype,error=error)
          DO i=1,50
             ii(i)=-1
          END DO
          stat=-1
          begin=1
          IF(ntype/=0) THEN
             DO
                IF (begin.eq.1) THEN
                   READ(parser%input_line,iostat=stat,FMT=*) itemp,itemp,ctemp,itemp,ftemp,ftemp,&
                        itemp,itemp,(ii(i),i=begin,ntype)
                ELSE IF (begin.gt.1) THEN
                   CALL parser_get_next_line(parser,1,error=error)
                   READ(parser%input_line,iostat=stat,FMT=*) (ii(i),i=begin,ntype)
                END IF
                DO i=ntype,1,-1
                   IF (ii(i).lt.0) THEN
                      begin=i
                   END IF
                END DO
                IF (stat.eq.0) EXIT
             END DO
          END IF
          ! 1-4 list
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,ntype,error=error)
          IF(ntype/=0) THEN
             itemp=(itemp-1)/6+1
             offset = 0
             IF(ASSOCIATED(conn_info%onfo_a)) offset = SIZE(conn_info%onfo_a)
             CALL reallocate(conn_info%onfo_a,1,offset+ntype)
             CALL reallocate(conn_info%onfo_b,1,offset+ntype)
             conn_info%onfo_a(offset+1:offset+ntype)=index_now
             DO i=1,50
                ii(i)=-1
             END DO
             stat=-1
             begin=1
             IF(ntype/=0) THEN
                DO
                   IF (begin.eq.1) THEN
                      READ(parser%input_line,iostat=stat,FMT=*) itemp,(ii(i),i=begin,ntype)
                   ELSE IF (begin.gt.1) THEN
                      CALL parser_get_next_line(parser,1,error=error)
                      READ(parser%input_line,iostat=stat,FMT=*) (ii(i),i=begin,ntype)
                   END IF
                   DO i=ntype,1,-1
                      IF (ii(i).lt.0) THEN
                         begin=i
                      END IF
                   END DO
                   IF (stat.eq.0) EXIT
                END DO
                DO i=1,ntype
                   conn_info%onfo_b(offset+i)=ii(i)
                END DO
             END IF
             CALL parser_get_next_line(parser,1,error=error)
          ELSE
             CALL parser_get_next_line(parser,1,error=error)
          END IF
       END DO
    END IF
    nsolute = natom
    CALL parser_release(parser,error=error)

    CALL parser_create(parser,file_name,para_env=para_env,error=error)
    ! BONDH SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the BONDH section'
    label = TRIM(avail_section(8))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ntype,error=error)
       offset = 0
       IF(ASSOCIATED(conn_info%bond_a)) offset = SIZE(conn_info%bond_a)
       CALL reallocate(conn_info%bond_a,1,offset+ntype)
       CALL reallocate(conn_info%bond_b,1,offset+ntype)
       CALL reallocate(conn_info%bond_type,1,offset+ntype)
       DO itype=1,ntype
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,conn_info%bond_a(offset+itype),error=error)
          CALL parser_get_object(parser,conn_info%bond_b(offset+itype),error=error)
          CALL parser_get_object(parser,itemp,error=error)
          conn_info%bond_type(offset+itype)=itemp
          IF(iw>0) WRITE(iw,*) "GTOP_INFO| PUT BONDH INFO HERE!!!!"
       END DO
       conn_info%bond_a(offset+1:offset+ntype) =  conn_info%bond_a(offset+1:offset+ntype) + natom_prev
       conn_info%bond_b(offset+1:offset+ntype) =  conn_info%bond_b(offset+1:offset+ntype) + natom_prev
    END IF
    ! BOND SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the BOND section'
    label = TRIM(avail_section(9))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ntype,error=error)
       offset = 0
       IF(ASSOCIATED(conn_info%bond_a)) offset = SIZE(conn_info%bond_a)
       CALL reallocate(conn_info%bond_a,1,offset+ntype)
       CALL reallocate(conn_info%bond_b,1,offset+ntype)
       CALL reallocate(conn_info%bond_type,1,offset+ntype)
       DO itype=1,ntype
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,conn_info%bond_a(offset+itype),error=error)
          CALL parser_get_object(parser,conn_info%bond_b(offset+itype),error=error)
          CALL parser_get_object(parser,itemp,error=error)
          conn_info%bond_type(offset+itype)=itemp
          IF(iw>0) WRITE(iw,*) "GTOP_INFO| PUT BOND INFO HERE!!!!"
       END DO
       conn_info%bond_a(offset+1:offset+ntype) =  conn_info%bond_a(offset+1:offset+ntype) + natom_prev
       conn_info%bond_b(offset+1:offset+ntype) =  conn_info%bond_b(offset+1:offset+ntype) + natom_prev
    END IF
    ! BONDANGLEH SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the BONDANGLEH section'
    label = TRIM(avail_section(11))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ntype,error=error)
       offset = 0
       IF(ASSOCIATED(conn_info%theta_a)) offset = SIZE(conn_info%theta_a)
       CALL reallocate(conn_info%theta_a,1,offset+ntype)
       CALL reallocate(conn_info%theta_b,1,offset+ntype)
       CALL reallocate(conn_info%theta_c,1,offset+ntype)
       CALL reallocate(conn_info%theta_type,1,offset+ntype)
       DO itype=1,ntype
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,conn_info%theta_a(offset+itype),error=error)
          CALL parser_get_object(parser,conn_info%theta_b(offset+itype),error=error)
          CALL parser_get_object(parser,conn_info%theta_c(offset+itype),error=error)
          CALL parser_get_object(parser,itemp,error=error)
          conn_info%theta_type(offset+itype)=itemp
          IF(iw>0) WRITE(iw,*) "GTOP_INFO| PUT BONDANGLEH INFO HERE!!!!"
       END DO
       conn_info%theta_a(offset+1:offset+ntype) =  conn_info%theta_a(offset+1:offset+ntype) + natom_prev
       conn_info%theta_b(offset+1:offset+ntype) =  conn_info%theta_b(offset+1:offset+ntype) + natom_prev
       conn_info%theta_c(offset+1:offset+ntype) =  conn_info%theta_c(offset+1:offset+ntype) + natom_prev
    END IF
    ! BONDANGLE SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the BONDANGLE section'
    label = TRIM(avail_section(12))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ntype,error=error)
       offset = 0
       IF(ASSOCIATED(conn_info%theta_a)) offset = SIZE(conn_info%theta_a)
       CALL reallocate(conn_info%theta_a,1,offset+ntype)
       CALL reallocate(conn_info%theta_b,1,offset+ntype)
       CALL reallocate(conn_info%theta_c,1,offset+ntype)
       CALL reallocate(conn_info%theta_type,1,offset+ntype)
       DO itype=1,ntype
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,conn_info%theta_a(offset+itype),error=error)
          CALL parser_get_object(parser,conn_info%theta_b(offset+itype),error=error)
          CALL parser_get_object(parser,conn_info%theta_c(offset+itype),error=error)
          CALL parser_get_object(parser,itemp,error=error)
          conn_info%theta_type(offset+itype)=itemp
          IF(iw>0) WRITE(iw,*) "GTOP_INFO| PUT BONDANGLE INFO HERE!!!!"
       END DO
       conn_info%theta_a(offset+1:offset+ntype) =  conn_info%theta_a(offset+1:offset+ntype) + natom_prev
       conn_info%theta_b(offset+1:offset+ntype) =  conn_info%theta_b(offset+1:offset+ntype) + natom_prev
       conn_info%theta_c(offset+1:offset+ntype) =  conn_info%theta_c(offset+1:offset+ntype) + natom_prev
    END IF
    ! IMPDIHEDRALH SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the IMPDIHEDRALH section'
    label = TRIM(avail_section(14))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ntype,error=error)
       offset = 0
       IF(ASSOCIATED(conn_info%impr_a)) offset = SIZE(conn_info%impr_a)
       CALL reallocate(conn_info%impr_a,1,offset+ntype)
       CALL reallocate(conn_info%impr_b,1,offset+ntype)
       CALL reallocate(conn_info%impr_c,1,offset+ntype)
       CALL reallocate(conn_info%impr_d,1,offset+ntype)
       CALL reallocate(conn_info%impr_type,1,offset+ntype)
       DO itype=1,ntype
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,conn_info%impr_a(offset+itype),error=error)
          CALL parser_get_object(parser,conn_info%impr_b(offset+itype),error=error)
          CALL parser_get_object(parser,conn_info%impr_c(offset+itype),error=error)
          CALL parser_get_object(parser,conn_info%impr_d(offset+itype),error=error)
          CALL parser_get_object(parser,itemp,error=error)
          conn_info%impr_type(offset+itype)=itemp
          IF(iw>0) WRITE(iw,*) "GTOP_INFO| PUT IMPDIHEDRALH INFO HERE!!!!"
       END DO
       conn_info%impr_a(offset+1:offset+ntype) =  conn_info%impr_a(offset+1:offset+ntype) + natom_prev
       conn_info%impr_b(offset+1:offset+ntype) =  conn_info%impr_b(offset+1:offset+ntype) + natom_prev
       conn_info%impr_c(offset+1:offset+ntype) =  conn_info%impr_c(offset+1:offset+ntype) + natom_prev
       conn_info%impr_d(offset+1:offset+ntype) =  conn_info%impr_d(offset+1:offset+ntype) + natom_prev
    END IF
    ! IMPDIHEDRAL SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the IMPDIHEDRAL section'
    label = TRIM(avail_section(15))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ntype,error=error)
       offset = 0
       IF(ASSOCIATED(conn_info%impr_a)) offset = SIZE(conn_info%impr_a)
       CALL reallocate(conn_info%impr_a,1,offset+ntype)
       CALL reallocate(conn_info%impr_b,1,offset+ntype)
       CALL reallocate(conn_info%impr_c,1,offset+ntype)
       CALL reallocate(conn_info%impr_d,1,offset+ntype)
       CALL reallocate(conn_info%impr_type,1,offset+ntype)
       DO itype=1,ntype
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,conn_info%impr_a(offset+itype),error=error)
          CALL parser_get_object(parser,conn_info%impr_b(offset+itype),error=error)
          CALL parser_get_object(parser,conn_info%impr_c(offset+itype),error=error)
          CALL parser_get_object(parser,conn_info%impr_d(offset+itype),error=error)
          CALL parser_get_object(parser,itemp,error=error)
          conn_info%impr_type(offset+itype)=itemp
          IF(iw>0) WRITE(iw,*) "GTOP_INFO| PUT IMPDIHEDRAL INFO HERE!!!!"
       END DO
       conn_info%impr_a(offset+1:offset+ntype) =  conn_info%impr_a(offset+1:offset+ntype) + natom_prev
       conn_info%impr_b(offset+1:offset+ntype) =  conn_info%impr_b(offset+1:offset+ntype) + natom_prev
       conn_info%impr_c(offset+1:offset+ntype) =  conn_info%impr_c(offset+1:offset+ntype) + natom_prev
       conn_info%impr_d(offset+1:offset+ntype) =  conn_info%impr_d(offset+1:offset+ntype) + natom_prev
    END IF
    ! DIHEDRALH SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the DIHEDRALH section'
    label = TRIM(avail_section(17))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ntype,error=error)
       offset = 0
       IF(ASSOCIATED(conn_info%phi_a)) offset = SIZE(conn_info%phi_a)
       CALL reallocate(conn_info%phi_a,1,offset+ntype)
       CALL reallocate(conn_info%phi_b,1,offset+ntype)
       CALL reallocate(conn_info%phi_c,1,offset+ntype)
       CALL reallocate(conn_info%phi_d,1,offset+ntype)
       CALL reallocate(conn_info%phi_type,1,offset+ntype)
       DO itype=1,ntype
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,conn_info%phi_a(offset+itype),error=error)
          CALL parser_get_object(parser,conn_info%phi_b(offset+itype),error=error)
          CALL parser_get_object(parser,conn_info%phi_c(offset+itype),error=error)
          CALL parser_get_object(parser,conn_info%phi_d(offset+itype),error=error)
          CALL parser_get_object(parser,itemp,error=error)
          conn_info%phi_type(offset+itype)=itemp
          IF(iw>0) WRITE(iw,*) "GTOP_INFO| PUT DIHEDRALH INFO HERE!!!!"
       END DO
       conn_info%phi_a(offset+1:offset+ntype) =  conn_info%phi_a(offset+1:offset+ntype) + natom_prev
       conn_info%phi_b(offset+1:offset+ntype) =  conn_info%phi_b(offset+1:offset+ntype) + natom_prev
       conn_info%phi_c(offset+1:offset+ntype) =  conn_info%phi_c(offset+1:offset+ntype) + natom_prev
       conn_info%phi_d(offset+1:offset+ntype) =  conn_info%phi_d(offset+1:offset+ntype) + natom_prev
    END IF
    ! DIHEDRAL SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the DIHEDRAL section'
    label = TRIM(avail_section(18))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ntype,error=error)
       offset = 0
       IF(ASSOCIATED(conn_info%phi_a)) offset = SIZE(conn_info%phi_a)
       CALL reallocate(conn_info%phi_a,1,offset+ntype)
       CALL reallocate(conn_info%phi_b,1,offset+ntype)
       CALL reallocate(conn_info%phi_c,1,offset+ntype)
       CALL reallocate(conn_info%phi_d,1,offset+ntype)
       CALL reallocate(conn_info%phi_type,1,offset+ntype)
       DO itype=1,ntype
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,conn_info%phi_a(offset+itype),error=error)
          CALL parser_get_object(parser,conn_info%phi_b(offset+itype),error=error)
          CALL parser_get_object(parser,conn_info%phi_c(offset+itype),error=error)
          CALL parser_get_object(parser,conn_info%phi_d(offset+itype),error=error)
          CALL parser_get_object(parser,itemp,error=error)
          conn_info%phi_type(offset+itype)=itemp
          IF(iw>0) WRITE(iw,*) "GTOP_INFO| PUT DIHEDRAL INFO HERE!!!!"
       END DO
       conn_info%phi_a(offset+1:offset+ntype) =  conn_info%phi_a(offset+1:offset+ntype) + natom_prev
       conn_info%phi_b(offset+1:offset+ntype) =  conn_info%phi_b(offset+1:offset+ntype) + natom_prev
       conn_info%phi_c(offset+1:offset+ntype) =  conn_info%phi_c(offset+1:offset+ntype) + natom_prev
       conn_info%phi_d(offset+1:offset+ntype) =  conn_info%phi_d(offset+1:offset+ntype) + natom_prev
    END IF
    CALL parser_release(parser,error=error)

    ! SOLVENTATOM and SOLVENTCONSTR SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the SOLVENTATOM section'
    nsolvent=( SIZE(atom_info%r(1,:))-nsolute ) / 3

    NULLIFY(na,am,ac,ba,bb)
    CALL parser_create(parser,file_name,para_env=para_env,error=error)
    label = TRIM(avail_section(20))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,natom,error=error)
       CALL reallocate(na,1,natom)
       CALL reallocate(am,1,natom)
       CALL reallocate(ac,1,natom)
       DO iatom=1,natom
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,itemp,error=error)
          CALL parser_get_object(parser,string,error=error)
          CALL parser_get_object(parser,na(iatom),error=error)
          CALL parser_get_object(parser,am(iatom),error=error)
          CALL parser_get_object(parser,ac(iatom),error=error)
          IF(iw>0) WRITE(iw,*) "GTOP_INFO| PUT SOLVENTATOM INFO HERE!!!!"
       END DO
    END IF
    label = TRIM(avail_section(21))
    ncon  = 0
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ncon,error=error)
       CALL reallocate(ba,1,ncon)
       CALL reallocate(bb,1,ncon)
       DO icon=1,ncon
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,ba(icon),error=error)
          CALL parser_get_object(parser,bb(icon),error=error)
       END DO
    END IF
    CALL parser_release(parser,error=error)

    offset = 0
    IF(ASSOCIATED(atom_info%id_molname)) offset=SIZE(atom_info%id_molname)
    CALL reallocate(atom_info%id_molname,1,offset+nsolvent*natom)
    CALL reallocate(atom_info%resid,1,offset+nsolvent*natom)
    CALL reallocate(atom_info%id_resname,1,offset+nsolvent*natom)
    CALL reallocate(atom_info%id_atmname,1,offset+nsolvent*natom)
    CALL reallocate(atom_info%id_element,1,offset+nsolvent*natom)
    CALL reallocate(atom_info%atm_charge,1,offset+nsolvent*natom)
    CALL reallocate(atom_info%atm_mass,1,offset+nsolvent*natom)
    DO isolvent=1,nsolvent
       offset = nsolute + natom*isolvent - natom
       DO iatom=1,natom
          index_now = offset+iatom
          atom_info%id_atmname(index_now) = str2id(s2s(namearray1(na(iatom))))
          atom_info%id_element(index_now) = str2id(s2s(namearray1(na(iatom))))
          atom_info%id_molname(index_now) = str2id(s2s("SOL"))
          atom_info%id_resname(index_now) = str2id(s2s("SOL"))
          atom_info%resid(index_now)   = isolvent
          atom_info%atm_mass(index_now)      = am(iatom)
          atom_info%atm_charge(index_now)    = ac(iatom)
       END DO
    END DO

    offset = 0
    IF(ASSOCIATED(conn_info%bond_a)) offset = SIZE(conn_info%bond_a)
    offset2= MAXVAL(conn_info%bond_type(:))
    CALL reallocate(conn_info%bond_a,1,offset+ncon*nsolvent)
    CALL reallocate(conn_info%bond_b,1,offset+ncon*nsolvent)
    CALL reallocate(conn_info%bond_type,1,offset+ncon*nsolvent)
    offset = offset - ncon
    DO isolvent=1,nsolvent
       offset = offset + ncon
       DO icon=1,ncon
          conn_info%bond_a(offset+icon)=nsolute+isolvent*ncon-ncon+ba(icon)
          conn_info%bond_b(offset+icon)=nsolute+isolvent*ncon-ncon+bb(icon)
          conn_info%bond_type(offset+icon)=offset2+isolvent*ncon-ncon+icon
       END DO
    END DO
    ! PARA_RES structure
    i = 0
    nbond_prev = 0
    IF(ASSOCIATED(conn_info%c_bond_a)) i = SIZE(conn_info%c_bond_a)
    nbond = SIZE(conn_info%bond_a)
    DO ibond=1+nbond_prev,nbond+nbond_prev
       iatom = conn_info%bond_a(ibond)
       jatom = conn_info%bond_b(ibond)
       IF(topology%para_res) THEN
          IF((atom_info%id_molname(iatom)/=atom_info%id_molname(jatom)).OR.&
               (atom_info%resid(iatom)/=atom_info%resid(jatom)).OR.&
               (atom_info%id_resname(iatom)/=atom_info%id_resname(jatom)))THEN
             IF(iw>0) WRITE(iw,'(T2,A,2I3)') "GTOP_INFO| PARA_RES, bond between molecules atom ",&
                  iatom,jatom
             i = i + 1
             CALL reallocate(conn_info%c_bond_a,1,i)
             CALL reallocate(conn_info%c_bond_b,1,i)
             CALL reallocate(conn_info%c_bond_type,1,i)
             conn_info%c_bond_a(i) = iatom
             conn_info%c_bond_b(i) = jatom
             conn_info%c_bond_type(i) = conn_info%bond_type(ibond)
          END IF
       ELSE
          IF(atom_info%id_molname(iatom)/=atom_info%id_molname(jatom)) THEN
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END IF
       END IF
    END DO

    DEALLOCATE(namearray1,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(namearray2,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(na,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(am,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(ac,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(ASSOCIATED(ba)) &
         DEALLOCATE(ba,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(ASSOCIATED(bb)) &
         DEALLOCATE(bb,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/GTOP_INFO",error=error)

  END SUBROUTINE read_topology_gromos

! *****************************************************************************
!> \brief ...
!> \param topology ...
!> \param para_env ...
!> \param subsys_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE read_coordinate_g96 (topology,para_env,subsys_section,error)
    TYPE(topology_parameters_type)           :: topology
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_coordinate_g96', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock = 1000

    CHARACTER(LEN=default_string_length)     :: label, string, strtmp, strtmp2
    CHARACTER(LEN=default_string_length), &
      DIMENSION(5)                           :: avail_section
    INTEGER                                  :: handle, istat, itemp, iw, &
                                                natom, newsize
    LOGICAL                                  :: failure, found
    REAL(KIND=dp)                            :: pfactor
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: velocity
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(section_vals_type), POINTER         :: velocity_section

    NULLIFY(parser, logger, velocity)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/G96_INFO",&
         extension=".subsysLog",error=error)
    CALL timeset(routineN,handle)

    pfactor=section_get_rval(subsys_section,"TOPOLOGY%MEMORY_PROGRESSION_FACTOR",error)
    atom_info => topology%atom_info
    IF (iw>0)  WRITE(iw,*) "    Reading in G96 file ",TRIM(topology%coord_file_name)
    avail_section( 1) = "TITLE"
    avail_section( 2) = "TIMESTEP"
    avail_section( 3) = "POSITION"
    avail_section( 4) = "VELOCITY"
    avail_section( 5) = "BOX"

    ! Element is assigned on the basis of the atm_name
    topology%aa_element = .TRUE.

    CALL reallocate(atom_info%id_molname,1,nblock)
    CALL reallocate(atom_info%id_resname,1,nblock)
    CALL reallocate(atom_info%resid,1,nblock)
    CALL reallocate(atom_info%id_atmname,1,nblock)
    CALL reallocate(atom_info%id_element,1,nblock)
    CALL reallocate(atom_info%r,1,3,1,nblock)
    CALL reallocate(atom_info%atm_mass,1,nblock)
    CALL reallocate(atom_info%atm_charge,1,nblock)
    CALL reallocate(atom_info%occup,1,nblock)
    CALL reallocate(atom_info%beta,1,nblock)
    ! TITLE SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'G96_INFO| Parsing the TITLE section'
    CALL parser_create(parser,topology%coord_file_name,para_env=para_env,error=error)
    label = TRIM(avail_section(1))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       DO
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,string,string_length=default_string_length,error=error)
          IF(string == TRIM("END")) EXIT
          IF(iw>0) WRITE(iw,*) "G96_INFO|   ",TRIM(string)
       END DO
    END IF
    CALL parser_release(parser,error=error)
    ! POSITION SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'G96_INFO| Parsing the POSITION section'
    CALL parser_create(parser,topology%coord_file_name,para_env=para_env,error=error)
    label = TRIM(avail_section(3))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       natom=0
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,string,string_length=default_string_length,error=error)
       DO
          IF (string == TRIM("END")) EXIT
          natom = natom + 1
          IF (natom > SIZE(atom_info%id_molname)) THEN
             newsize = INT(pfactor*natom)
             CALL reallocate(atom_info%id_molname,1,newsize)
             CALL reallocate(atom_info%id_resname,1,newsize)
             CALL reallocate(atom_info%resid,1,newsize)
             CALL reallocate(atom_info%id_atmname,1,newsize)
             CALL reallocate(atom_info%id_element,1,newsize)
             CALL reallocate(atom_info%r,1,3,1,newsize)
             CALL reallocate(atom_info%atm_mass,1,newsize)
             CALL reallocate(atom_info%atm_charge,1,newsize)
             CALL reallocate(atom_info%occup,1,newsize)
             CALL reallocate(atom_info%beta,1,newsize)
          END IF
          READ (string,*) &
               atom_info%resid(natom),strtmp,strtmp2,&
               itemp, atom_info%r(1,natom), atom_info%r(2,natom), atom_info%r(3,natom)
          atom_info%id_resname(natom)=str2id(s2s(strtmp))
          atom_info%id_atmname(natom)=str2id(s2s(strtmp2))
          atom_info%id_molname(natom) = atom_info%id_resname(natom)
          atom_info%id_element(natom) = atom_info%id_atmname(natom)
          atom_info%r(1,natom) = cp_unit_to_cp2k(atom_info%r(1,natom),"nm",error=error)
          atom_info%r(2,natom) = cp_unit_to_cp2k(atom_info%r(2,natom),"nm",error=error)
          atom_info%r(3,natom) = cp_unit_to_cp2k(atom_info%r(3,natom),"nm",error=error)
          IF(iw>0) WRITE(iw,*) "G96_INFO| PUT POSITION INFO HERE!!!!"
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,string,string_length=default_string_length,error=error)
       END DO
    END IF
    CALL parser_release(parser,error=error)
    CALL reallocate(velocity,1,3,1,natom)
    ! VELOCITY SECTION
    IF(topology%use_g96_velocity) THEN
       IF(iw>0) WRITE(iw,'(T2,A)') 'G96_INFO| Parsing the VELOCITY section'
       CALL parser_create(parser,topology%coord_file_name,para_env=para_env,error=error)
       label = TRIM(avail_section(4))
       CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
       IF(found) THEN
          natom=0
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,string,string_length=default_string_length,error=error)
          DO
             IF(string == TRIM("END")) EXIT
             natom=natom+1
             READ (string,*) &
                  atom_info%resid(natom),strtmp,strtmp2,&
                  itemp,velocity(1,natom),velocity(2,natom),velocity(3,natom)
             atom_info%id_resname(natom)=str2id(strtmp)
             atom_info%id_atmname(natom)=str2id(strtmp2)
             atom_info%id_molname(natom)=atom_info%id_resname(natom)
             atom_info%id_element(natom)=atom_info%id_atmname(natom)
             velocity(1,natom) = cp_unit_to_cp2k(velocity(1,natom),"nm*ps^-1",error=error)
             velocity(2,natom) = cp_unit_to_cp2k(velocity(2,natom),"nm*ps^-1",error=error)
             velocity(3,natom) = cp_unit_to_cp2k(velocity(3,natom),"nm*ps^-1",error=error)
             IF(iw>0) WRITE(iw,*) "G96_INFO| PUT VELOCITY INFO HERE!!!!"
             CALL parser_get_next_line(parser,1,error=error)
             CALL parser_get_object(parser,string,string_length=default_string_length,error=error)
          END DO
          CALL parser_release(parser,error=error)
          velocity_section => section_vals_get_subs_vals(subsys_section,"VELOCITY",error=error)
          CALL section_velocity_val_set(velocity_section, velocity=velocity, &
               conv_factor=1.0_dp, error=error)
       END IF
    END IF
    DEALLOCATE(velocity,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL reallocate(atom_info%id_molname,1,natom)
    CALL reallocate(atom_info%id_resname,1,natom)
    CALL reallocate(atom_info%resid,1,natom)
    CALL reallocate(atom_info%id_atmname,1,natom)
    CALL reallocate(atom_info%id_element,1,natom)
    CALL reallocate(atom_info%r,1,3,1,natom)
    CALL reallocate(atom_info%atm_mass,1,natom)
    CALL reallocate(atom_info%atm_charge,1,natom)
    CALL reallocate(atom_info%occup,1,natom)
    CALL reallocate(atom_info%beta,1,natom)

    IF(.NOT.topology%para_res) atom_info%resid(:) = 1

    topology%natoms = natom
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/G96_INFO",error=error)
    CALL timestop(handle)

  END SUBROUTINE read_coordinate_g96

END MODULE topology_gromos

