!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
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
!!     IKUO 08.01.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE topology_psf
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp,&
                                             default_string_length
  USE memory_utilities,                ONLY: reallocate
  USE cp_parser,                       ONLY: parser_release,&
                                             parser_get_object,&
                                             parser_get_next_line,&
                                             parser_search_string,&
                                             parser_create,&
                                             cp_parser_type
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: topology_parameters_type,&
                                             atom_info_type,&
                                             connectivity_info_type
  USE string_utilities,                ONLY: uppercase
  USE cp_output_handling,              ONLY: cp_print_key_unit_nr, &
                                             cp_print_key_finished_output  

#include "cp_common_uses.h"
  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "topology_psf"

  PRIVATE
  PUBLIC :: read_topology_psf, write_topology_psf

!!*****
!******************************************************************************

CONTAINS

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

  SUBROUTINE read_topology_psf (topology,globenv, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      pointer, OPTIONAL                      :: globenv
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error 
    ! Local Variables
    TYPE (atom_info_type),POINTER            :: atom_info
    TYPE (connectivity_info_type),POINTER    :: conn_info

    CHARACTER(LEN=default_string_length)              :: label, dummy_field
    INTEGER                                  :: i, iatom, ibond, stat, &
                                                ilen, iphi, itheta, iunit, &
                                                iw, natom, nbond, nphi, ntheta
    CHARACTER(len=*), PARAMETER :: routineN = 'read_coordinate_psf', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: handle
    LOGICAL                                  :: found, failure
    INTEGER :: jatom
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(parser, logger)
    logger => cp_error_get_logger(error)
    CALL write_checkpoint_information("Entering "//routineN,globenv,error)
    iw = cp_print_key_unit_nr(logger,globenv%input_file,"FORCE_EVAL%MM%PRINT%TOPOLOGY_INFO/PSF_INFO",&
         extension=".mmLog",error=error)    
    CALL timeset(routineN,'I','',handle)
    CALL parser_create(parser,topology%conn_file_name,para_env=globenv%para_env)
    
    atom_info => topology%atom_info
    conn_info => topology%conn_info
    
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
       IF (ASSOCIATED(atom_info%label_molname)) DEALLOCATE(atom_info%label_molname,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(atom_info%label_molname(natom),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     
       IF (ASSOCIATED(atom_info%label_resid)) DEALLOCATE(atom_info%label_resid,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(atom_info%label_resid(natom),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     
       IF (ASSOCIATED(atom_info%label_resname)) DEALLOCATE(atom_info%label_resname,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(atom_info%label_resname(natom),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       
       IF (ASSOCIATED(atom_info%label_atmname)) DEALLOCATE(atom_info%label_atmname,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(atom_info%label_atmname(natom),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     
       IF (ASSOCIATED(atom_info%atm_charge)) DEALLOCATE(atom_info%atm_charge,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(atom_info%atm_charge(natom),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     
       IF (ASSOCIATED(atom_info%atm_mass)) DEALLOCATE(atom_info%atm_mass,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(atom_info%atm_mass(natom),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     
       !Read in the atom info
       DO iatom=1,natom
          CALL parser_get_object(parser,i,newline=.TRUE.)
          CALL parser_get_object(parser,atom_info%label_molname(iatom))
          CALL parser_get_object(parser,atom_info%label_resid(iatom))
          CALL parser_get_object(parser,atom_info%label_resname(iatom))
          CALL parser_get_object(parser,dummy_field)
          CALL parser_get_object(parser,atom_info%label_atmname(iatom))
          CALL parser_get_object(parser,atom_info%atm_charge(iatom))
          CALL parser_get_object(parser,atom_info%atm_mass(iatom))
          CALL uppercase ( atom_info%label_atmname(iatom) )
       END DO
    END IF

    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
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
       ALLOCATE(conn_info%bond_a(nbond),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(conn_info%bond_b(nbond),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       !Read in the atom info
       DO ibond=1,nbond,4
          CALL parser_get_next_line(parser,1,error=error)
          i=0
          DO WHILE ((i<4).AND.((ibond+i)<=nbond))
             CALL parser_get_object(parser,conn_info%bond_a(ibond+i))
             CALL parser_get_object(parser,conn_info%bond_b(ibond+i))
             i=i+1
          END DO
       END DO
    END IF
    
    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    i = 0
    DO ibond=1,nbond
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

    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
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
       ALLOCATE(conn_info%theta_a(ntheta),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(conn_info%theta_b(ntheta),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(conn_info%theta_c(ntheta),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       !Read in the atom info
       DO itheta=1,ntheta,3
          CALL parser_get_next_line(parser,1,error=error)
          i=0
          DO WHILE ((i<3).AND.((itheta+i)<=ntheta))
             CALL parser_get_object(parser,conn_info%theta_a(itheta+i))
             CALL parser_get_object(parser,conn_info%theta_b(itheta+i))
             CALL parser_get_object(parser,conn_info%theta_c(itheta+i))
             i=i+1
          END DO
       END DO
    END IF

    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    ALLOCATE(conn_info%ub_a(ntheta),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(conn_info%ub_b(ntheta),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(conn_info%ub_c(ntheta),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    conn_info%ub_a(:) = conn_info%theta_a(:)
    conn_info%ub_b(:) = conn_info%theta_b(:)
    conn_info%ub_c(:) = conn_info%theta_c(:)

    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    IF(iw>0) WRITE(iw,*) 'Parsing the NPHI (ie dihedrals) section'
    label = '!NPHI'
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF ( .NOT. found ) THEN
       IF(iw>0) WRITE ( iw, '( A )' ) ' No NPHI section '
       nphi = 0
    ELSE
       CALL parser_get_object  (parser,nphi)
       IF(iw>0) WRITE(iw,*) 'NPHI = ',nphi
       !malloc the memory that we need
       ALLOCATE(conn_info%phi_a(nphi),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(conn_info%phi_b(nphi),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(conn_info%phi_c(nphi),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(conn_info%phi_d(nphi),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       !Read in the atom info
       DO iphi=1,nphi,2
          CALL parser_get_next_line(parser,1,error=error)
          i=0
          DO WHILE ((i<2).AND.((iphi+i)<=nphi))
             CALL parser_get_object(parser,conn_info%phi_a(iphi+i))
             CALL parser_get_object(parser,conn_info%phi_b(iphi+i))
             CALL parser_get_object(parser,conn_info%phi_c(iphi+i))
             CALL parser_get_object(parser,conn_info%phi_d(iphi+i))
             i=i+1
          END DO
       END DO
    END IF

    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    ALLOCATE(conn_info%onfo_a(nphi),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(conn_info%onfo_b(nphi),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    conn_info%onfo_a(:) = conn_info%phi_a(:)
    conn_info%onfo_b(:) = conn_info%phi_d(:)

    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    IF(iw>0) WRITE(iw,*) 'Parsing the NIMPHI (ie dihedrals) section'
    label = '!NIMPHI'
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF ( .NOT. found ) THEN
       IF(iw>0) WRITE ( iw, '( A )' ) ' No NIMPHI section '
       nphi = 0
    ELSE
       CALL parser_get_object  (parser,nphi)
       IF(iw>0) WRITE(iw,*) 'NIMPR = ',nphi
       !malloc the memory that we need
       ALLOCATE(conn_info%impr_a(nphi),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(conn_info%impr_b(nphi),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)       
       ALLOCATE(conn_info%impr_c(nphi),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(conn_info%impr_d(nphi),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       !Read in the atom info
       DO iphi=1,nphi,2
          CALL parser_get_next_line(parser,1,error=error)
          i=0
          DO WHILE ((i<2).AND.((iphi+i)<=nphi))
             CALL parser_get_object(parser,conn_info%impr_a(iphi+i))
             CALL parser_get_object(parser,conn_info%impr_b(iphi+i))
             CALL parser_get_object(parser,conn_info%impr_c(iphi+i))
             CALL parser_get_object(parser,conn_info%impr_d(iphi+i))
             i=i+1
          END DO
       END DO
    END IF

    CALL parser_release(parser,error=error)
    CALL timestop(0.0_dp,handle)
    CALL write_checkpoint_information("exiting "//routineN,globenv,error)
    CALL cp_print_key_finished_output(iw,logger,globenv%input_file,&
         "FORCE_EVAL%MM%PRINT%TOPOLOGY_INFO/PSF_INFO",error=error)

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

  SUBROUTINE write_topology_psf (topology,globenv,error)
    IMPLICIT NONE
    !Arguments
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      pointer, OPTIONAL                   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error 
    ! Local Variables
    INTEGER                                  :: i, j, iatom, ibond, &
                                                ilen, iphi, itheta, iunit, &
                                                iw, natom, nbond, nphi, ntheta
    CHARACTER(len=*), PARAMETER :: routineN = 'write_topology_psf', &
      routineP = moduleN//':'//routineN
    TYPE (atom_info_type),POINTER            :: atom_info
    TYPE (connectivity_info_type),POINTER    :: conn_info
    CHARACTER(LEN=default_string_length)     :: record
    INTEGER                                  :: file_unit,imol,id1,id2,handle
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    CALL write_checkpoint_information("Entering "//routineN,globenv,error)
    iw = cp_print_key_unit_nr(logger,globenv%input_file,"FORCE_EVAL%MM%PRINT%TOPOLOGY_INFO/PSF_INFO",&
         extension=".mmLog",error=error)    
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
       WRITE(file_unit,'(I8,1X,A5,I4,1X,A5,A5,A5,F10.6,4X,F10.4,4X,I8)') &
            i,&
            atom_info%label_molname(i),&
            atom_info%map_mol_num(i),&
            atom_info%label_resname(i),&
            atom_info%label_atmname(i),&
            atom_info%label_atmname(i),&
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
    CALL write_checkpoint_information("exiting "//routineN,globenv,error)
    CALL cp_print_key_finished_output(iw,logger,globenv%input_file,&
         "FORCE_EVAL%MM%PRINT%TOPOLOGY_INFO/PSF_INFO",error=error)  
    
  END SUBROUTINE write_topology_psf

END MODULE topology_psf

!******************************************************************************
