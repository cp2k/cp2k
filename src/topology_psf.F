!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****** cp2k/topology_psf [1.0] *
!!
!!   NAME
!!     topology_psf
!!
!!   FUNCTION
!!     Functionality to read in PSF topologies and convert it into local
!!     data structures
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE topology_psf
  USE cp_files,                        ONLY: open_file
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
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_psf'

  PRIVATE
  PUBLIC :: pre_read_topology_psf, read_topology_psf, write_topology_psf

!!*****
!******************************************************************************

CONTAINS

!!*****
!******************************************************************************
!!****** topology/pre_read_topology_psf [1.0] *
!!
!!   NAME
!!     pre_read_topology_psf
!!
!!   FUNCTION
!!
!!   SOURCE
!******************************************************************************

  SUBROUTINE pre_read_topology_psf (topology,subsys_section, error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pre_read_topology_psf', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw, stat
    LOGICAL                                  :: failure
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/PSF_INFO",&
         extension=".subsysLog",error=error)    
    CALL timeset(routineN,'I','',handle)
    
    atom_info => topology%atom_info
    conn_info => topology%conn_info

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
         "PRINT%TOPOLOGY_INFO/PSF_INFO",error=error)

  END SUBROUTINE pre_read_topology_psf

!!*****
!******************************************************************************
!!****** topology/read_topology_psf [1.0] *
!!
!!   NAME
!!     read_topology_psf
!!
!!   FUNCTION
!!     Read PSF topology file
!!
!!   SOURCE
!******************************************************************************

  SUBROUTINE read_topology_psf (filename,topology,globenv,subsys_section, error)
    CHARACTER(LEN=default_string_length), &
      INTENT(IN)                             :: filename
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_topology_psf', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: dummy_field, label
    INTEGER :: handle, i, iatom, ibond, index_now, iphi, itheta, iw, jatom, &
      natom, natom_prev, nbond, nbond_prev, nphi, nphi_prev, ntheta, &
      ntheta_prev
    LOGICAL                                  :: failure, found
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_parser_type), POINTER            :: parser

    failure = .FALSE.
    NULLIFY(parser, logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/PSF_INFO",&
         extension=".subsysLog",error=error)    
    CALL timeset(routineN,'I','',handle)
    CALL parser_create(parser,filename,para_env=globenv%para_env)
    
    atom_info => topology%atom_info
    conn_info => topology%conn_info

    natom_prev = 0
    IF(ASSOCIATED(atom_info%label_molname)) natom_prev = SIZE(atom_info%label_molname)

    IF(iw>0) WRITE(iw,*) 'Parsing the NATOM section'

    label = '!NATOM'
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF ( .NOT. found ) THEN
       IF(iw>0) WRITE ( iw, '( A )' ) ' No NATOM section '
       natom = 0
    ELSE
       CALL parser_get_object  (parser,natom)
       IF(iw>0) WRITE(iw,*) 'NATOM = ',natom
       !malloc the memory that we need
       CALL reallocate(atom_info%label_molname,1,natom_prev+natom)
       CALL reallocate(atom_info%label_resid,1,natom_prev+natom)
       CALL reallocate(atom_info%label_resname,1,natom_prev+natom)
       CALL reallocate(atom_info%label_atmname,1,natom_prev+natom)
       CALL reallocate(atom_info%atm_charge,1,natom_prev+natom)
       CALL reallocate(atom_info%atm_mass,1,natom_prev+natom)
     
       !Read in the atom info
       DO iatom=1,natom
          index_now=iatom+natom_prev
          CALL parser_get_object(parser,i,newline=.TRUE.)
          CALL parser_get_object(parser,atom_info%label_molname(index_now))
          CALL parser_get_object(parser,atom_info%label_resid(index_now))
          CALL parser_get_object(parser,atom_info%label_resname(index_now))
          CALL parser_get_object(parser,dummy_field)
          CALL parser_get_object(parser,atom_info%label_atmname(index_now))
          CALL parser_get_object(parser,atom_info%atm_charge(index_now))
          CALL parser_get_object(parser,atom_info%atm_mass(index_now))
          CALL uppercase ( atom_info%label_atmname(index_now) )
       END DO
    END IF

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    nbond_prev = 0
    IF(ASSOCIATED(conn_info%bond_a)) nbond_prev = SIZE(conn_info%bond_a)

    IF(iw>0) WRITE(iw,*) 'Parsing the NBOND section'
    label = '!NBOND'
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF ( .NOT. found ) THEN
       IF(iw>0) WRITE ( iw, '( A )' ) ' No NBOND section '
       nbond = 0
    ELSE
       CALL parser_get_object  (parser,nbond)
       
       IF(iw>0) WRITE(iw,*) 'NBOND = ',nbond
       !malloc the memory that we need
       CALL reallocate(conn_info%bond_a,1,nbond_prev+nbond)
       CALL reallocate(conn_info%bond_b,1,nbond_prev+nbond)
       !Read in the atom info
       DO ibond=1,nbond,4
          CALL parser_get_next_line(parser,1,error=error)
          i=0
          DO WHILE ((i<4).AND.((ibond+i)<=nbond))
             index_now=nbond_prev+ibond+i
             CALL parser_get_object(parser,conn_info%bond_a(index_now))
             CALL parser_get_object(parser,conn_info%bond_b(index_now))
             conn_info%bond_a(index_now)=conn_info%bond_a(index_now)+natom_prev
             conn_info%bond_b(index_now)=conn_info%bond_b(index_now)+natom_prev
             i=i+1
          END DO
       END DO
    END IF
    
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    i = 0
    IF(ASSOCIATED(conn_info%c_bond_a)) i = SIZE(conn_info%c_bond_a)
    DO ibond=1+nbond_prev,nbond+nbond_prev
       iatom = conn_info%bond_a(ibond)
       jatom = conn_info%bond_b(ibond)
       IF(topology%para_res) THEN
          IF((atom_info%label_molname(iatom)/=atom_info%label_molname(jatom)).OR.&
               (atom_info%label_resid(iatom)/=atom_info%label_resid(jatom)).OR.&
               (atom_info%label_resname(iatom)/=atom_info%label_resname(jatom)))THEN
             IF(iw>0) WRITE(iw,*) "      PARA_RES, bond between molecules atom ",&
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

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ntheta_prev = 0
    IF(ASSOCIATED(conn_info%theta_a)) ntheta_prev = SIZE(conn_info%theta_a)

    IF(iw>0) WRITE(iw,*) 'Parsing the NTHETA section'
    label = '!NTHETA'
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF ( .NOT. found ) THEN
       IF(iw>0) WRITE ( iw, '( A )' ) ' No NTHETA section '
       ntheta = 0
    ELSE
       CALL parser_get_object  (parser,ntheta)
       IF(iw>0) WRITE(iw,*) 'NTHETA = ',ntheta
       !malloc the memory that we need
       CALL reallocate(conn_info%theta_a,1,ntheta_prev+ntheta)
       CALL reallocate(conn_info%theta_b,1,ntheta_prev+ntheta)
       CALL reallocate(conn_info%theta_c,1,ntheta_prev+ntheta)
       CALL reallocate(conn_info%ub_a,1,ntheta_prev+ntheta)
       CALL reallocate(conn_info%ub_b,1,ntheta_prev+ntheta)
       CALL reallocate(conn_info%ub_c,1,ntheta_prev+ntheta)
       !Read in the atom info
       DO itheta=1,ntheta,3
          CALL parser_get_next_line(parser,1,error=error)
          i=0
          DO WHILE ((i<3).AND.((itheta+i)<=ntheta))
             index_now=ntheta_prev+itheta+i
             CALL parser_get_object(parser,conn_info%theta_a(index_now))
             CALL parser_get_object(parser,conn_info%theta_b(index_now))
             CALL parser_get_object(parser,conn_info%theta_c(index_now))
             conn_info%theta_a(index_now)=conn_info%theta_a(index_now)+natom_prev
             conn_info%theta_b(index_now)=conn_info%theta_b(index_now)+natom_prev
             conn_info%theta_c(index_now)=conn_info%theta_c(index_now)+natom_prev
             i=i+1
          END DO
       END DO
    END IF

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    conn_info%ub_a(:) = conn_info%theta_a(:)
    conn_info%ub_b(:) = conn_info%theta_b(:)
    conn_info%ub_c(:) = conn_info%theta_c(:)
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    nphi_prev = 0
    IF(ASSOCIATED(conn_info%phi_a)) nphi_prev = SIZE(conn_info%phi_a)

    IF(iw>0) WRITE(iw,*) 'Parsing the NPHI section'
    label = '!NPHI'
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF ( .NOT. found ) THEN
       IF(iw>0) WRITE ( iw, '( A )' ) ' No NPHI section '
       nphi = 0
    ELSE
       CALL parser_get_object  (parser,nphi)
       IF(iw>0) WRITE(iw,*) 'NPHI = ',nphi
       !malloc the memory that we need
       CALL reallocate(conn_info%phi_a,1,nphi_prev+nphi)
       CALL reallocate(conn_info%phi_b,1,nphi_prev+nphi)
       CALL reallocate(conn_info%phi_c,1,nphi_prev+nphi)
       CALL reallocate(conn_info%phi_d,1,nphi_prev+nphi)
       CALL reallocate(conn_info%onfo_a,1,nphi_prev+nphi)
       CALL reallocate(conn_info%onfo_b,1,nphi_prev+nphi)
       !Read in the atom info
       DO iphi=1,nphi,2
          CALL parser_get_next_line(parser,1,error=error)
          i=0
          DO WHILE ((i<2).AND.((iphi+i)<=nphi))
             index_now=nphi_prev+iphi+i
             CALL parser_get_object(parser,conn_info%phi_a(index_now))
             CALL parser_get_object(parser,conn_info%phi_b(index_now))
             CALL parser_get_object(parser,conn_info%phi_c(index_now))
             CALL parser_get_object(parser,conn_info%phi_d(index_now))
             conn_info%phi_a(index_now)=conn_info%phi_a(index_now)+natom_prev
             conn_info%phi_b(index_now)=conn_info%phi_b(index_now)+natom_prev
             conn_info%phi_c(index_now)=conn_info%phi_c(index_now)+natom_prev
             conn_info%phi_d(index_now)=conn_info%phi_d(index_now)+natom_prev
             i=i+1
          END DO
       END DO
    END IF

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    conn_info%onfo_a(:) = conn_info%phi_a(:)
    conn_info%onfo_b(:) = conn_info%phi_d(:)
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    nphi_prev = 0
    IF(ASSOCIATED(conn_info%impr_a)) nphi_prev = SIZE(conn_info%impr_a)

    IF(iw>0) WRITE(iw,*) 'Parsing the NIMPHI section'
    label = '!NIMPHI'
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF ( .NOT. found ) THEN
       IF(iw>0) WRITE ( iw, '( A )' ) ' No NIMPHI section '
       nphi = 0
    ELSE
       CALL parser_get_object  (parser,nphi)
       IF(iw>0) WRITE(iw,*) 'NIMPR = ',nphi
       !malloc the memory that we need
       CALL reallocate(conn_info%impr_a,1,nphi_prev+nphi)
       CALL reallocate(conn_info%impr_b,1,nphi_prev+nphi)
       CALL reallocate(conn_info%impr_c,1,nphi_prev+nphi)
       CALL reallocate(conn_info%impr_d,1,nphi_prev+nphi)
       !Read in the atom info
       DO iphi=1,nphi,2
          CALL parser_get_next_line(parser,1,error=error)
          i=0
          DO WHILE ((i<2).AND.((iphi+i)<=nphi))
             index_now=nphi_prev+iphi+i
             CALL parser_get_object(parser,conn_info%impr_a(index_now))
             CALL parser_get_object(parser,conn_info%impr_b(index_now))
             CALL parser_get_object(parser,conn_info%impr_c(index_now))
             CALL parser_get_object(parser,conn_info%impr_d(index_now))
             conn_info%impr_a(index_now)=conn_info%impr_a(index_now)+natom_prev
             conn_info%impr_b(index_now)=conn_info%impr_b(index_now)+natom_prev
             conn_info%impr_c(index_now)=conn_info%impr_c(index_now)+natom_prev
             conn_info%impr_d(index_now)=conn_info%impr_d(index_now)+natom_prev
             i=i+1
          END DO
       END DO
    END IF

    CALL parser_release(parser,error=error)
    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/PSF_INFO",error=error)

  END SUBROUTINE read_topology_psf

!!*****
!******************************************************************************
!!****** topology/write_topology_psf [1.0] *
!!
!!   NAME
!!     write_topology_psf
!!
!!   FUNCTION
!!     Write PSF topology file
!!
!!   SOURCE
!******************************************************************************

  SUBROUTINE write_topology_psf (topology,globenv,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_topology_psf', &
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
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/PSF_INFO",&
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
       WRITE(file_unit,'(I8,1X,A5,I4,1X,A5,A5,A5,F10.6,4X,F10.4,4X,I8)') &
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

    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/PSF_INFO",error=error)  
    
  END SUBROUTINE write_topology_psf

END MODULE topology_psf

!******************************************************************************
