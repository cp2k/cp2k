!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Functionality to read in PSF topologies and convert it into local
!>      data structures
!> \author ikuo
!>      tlaino 10.2006
! *****************************************************************************
MODULE topology_psf
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_generate_filename,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                             parser_get_object,&
                                             parser_search_string,&
                                             parser_test_next_token
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE force_fields_input,              ONLY: read_chrg_section
  USE input_constants,                 ONLY: do_conn_psf,&
                                             do_conn_psf_u
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
  USE string_table,                    ONLY: id2str,&
                                             s2s,&
                                             str2id
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: atom_info_type,&
                                             connectivity_info_type,&
                                             topology_parameters_type
  USE topology_util,                   ONLY: array1_list_type,&
                                             reorder_structure,&
                                             tag_molecule
  USE util,                            ONLY: sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_psf'

  PRIVATE
  PUBLIC :: read_topology_psf,&
            write_topology_psf,&
            psf_post_process,&
            idm_psf

CONTAINS

! *****************************************************************************
!> \brief Read PSF topology file
!>      Teodoro Laino - Introduced CHARMM31 EXT PSF standard format
!> \param filename ...
!> \param topology ...
!> \param para_env ...
!> \param subsys_section ...
!> \param psf_type ...
!> \param error ...
!> \par History
!>      04-2007 Teodoro Laino - Zurich University [tlaino]
!>      This routine should contain only information read from the PSF format
!>      and all post_process should be performef in the psf_post_process
! *****************************************************************************
  SUBROUTINE read_topology_psf (filename,topology,para_env,subsys_section, psf_type, error)
    CHARACTER(LEN=*), INTENT(IN)             :: filename
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    INTEGER, INTENT(IN)                      :: psf_type
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_topology_psf', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2*default_string_length)   :: psf_format
    CHARACTER(LEN=3)                         :: c_int
    CHARACTER(LEN=default_string_length)     :: dummy_field, field, label, &
                                                strtmp1, strtmp2, strtmp3
    INTEGER :: handle, i, iatom, ibond, idum, index_now, iphi, itheta, iw, &
      natom, natom_prev, nbond, nbond_prev, nphi, nphi_prev, ntheta, &
      ntheta_prev, output_unit
    LOGICAL                                  :: failure, found
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_parser_type), POINTER            :: parser

    NULLIFY(parser, logger)
    failure = .FALSE.
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/PSF_INFO",&
         extension=".subsysLog",error=error)
    CALL timeset(routineN,handle)
    CALL parser_create(parser,filename,para_env=para_env,error=error)

    atom_info => topology%atom_info
    conn_info => topology%conn_info
    natom_prev = 0
    IF(ASSOCIATED(atom_info%id_molname)) natom_prev = SIZE(atom_info%id_molname)
    c_int      = 'I8'
    label      = 'PSF'
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF (.NOT. found) THEN
       IF (output_unit>0) THEN
          WRITE(output_unit,'(A)')"ERROR| Missing PSF specification line"
       END IF
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF
    DO WHILE (parser_test_next_token(parser,error=error)/="EOL")
       CALL parser_get_object  (parser,field,error=error)
       SELECT CASE(field(1:3))
       CASE ("PSF")
          IF (psf_type==do_conn_psf) THEN
            ! X-PLOR PSF format "similar" to the plain CHARMM PSF format
            psf_format = '(I8,1X,A4,I5,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8)'
          ENDIF
       CASE ("EXT")
          IF (psf_type==do_conn_psf) THEN
            ! EXTEnded CHARMM31 format
            psf_format = '(I10,T12,A7,T21,I8,T30,A7,T39,A6,T47,A6,T53,F10.6,T69,F8.3,T88,I1)'
            c_int      = 'I10'
          ELSE
            CALL stop_program(routineN,moduleN,__LINE__,&
                              "PSF_INFO| "//field(1:3)//" :: not available for UPSF format!")
          ENDIF
       CASE DEFAULT
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "PSF_INFO| "//field(1:3)//" :: Unimplemented keyword in CP2K PSF/UPSF format!")
       END SELECT
    END DO
    IF(iw>0) WRITE(iw,'(T2,A)') 'PSF_INFO| Parsing the NATOM section'
    !
    ! ATOM section
    !
    label = '!NATOM'
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF ( .NOT. found ) THEN
       IF(iw>0) WRITE (iw,'(T2,A)') 'PSF_INFO| No NATOM section '
       natom = 0
    ELSE
       CALL parser_get_object  (parser,natom,error=error)
       CALL cp_assert(natom_prev+natom<=topology%natoms,&
            cp_failure_level,cp_assertion_failed,routineP,&
            "Number of atoms in connectivity control is larger than the "//&
            "number of atoms in coordinate control. check coordinates and "//&
            "connectivity. "//&
CPSourceFileRef,&
             error=error,failure=failure)
       IF(iw>0) WRITE(iw,'(T2,A,'//TRIM(c_int)//')') 'PSF_INFO| NATOM = ',natom
       !malloc the memory that we need
       CALL reallocate(atom_info%id_molname,1,natom_prev+natom)
       CALL reallocate(atom_info%resid,1,natom_prev+natom)
       CALL reallocate(atom_info%id_resname,1,natom_prev+natom)
       CALL reallocate(atom_info%id_atmname,1,natom_prev+natom)
       CALL reallocate(atom_info%atm_charge,1,natom_prev+natom)
       CALL reallocate(atom_info%atm_mass,1,natom_prev+natom)
       !Read in the atom info
       IF (psf_type==do_conn_psf_u) THEN
          DO iatom=1,natom
             index_now=iatom+natom_prev
             CALL parser_get_next_line(parser,1,error=error)
             READ(parser%input_line,FMT=*,ERR=9)i,&
                                      strtmp1,&
                                      atom_info%resid(index_now),&
                                      strtmp2,&
                                      dummy_field,&
                                      strtmp3,&
                                      atom_info%atm_charge(index_now),&
                                      atom_info%atm_mass(index_now)
             atom_info%id_molname(index_now)=str2id(s2s(strtmp1))
             atom_info%id_resname(index_now)=str2id(s2s(strtmp2))
             atom_info%id_atmname(index_now)=str2id(s2s(strtmp3))
          END DO
       ELSE
          DO iatom=1,natom
             index_now=iatom+natom_prev
             CALL parser_get_next_line(parser,1,error=error)
             READ(parser%input_line,FMT=psf_format)&
                  i,&
                  strtmp1,&
                  atom_info%resid(index_now),&
                  strtmp2,&
                  dummy_field,&
                  strtmp3,&
                  atom_info%atm_charge(index_now),&
                  atom_info%atm_mass(index_now),&
                  idum
             atom_info%id_molname(index_now)=str2id(s2s(strtmp1))
             atom_info%id_resname(index_now)=str2id(s2s(strtmp2))
             atom_info%id_atmname(index_now)=str2id(s2s(ADJUSTL(strtmp3)))
          END DO
       END IF
    END IF

    !
    ! BOND section
    !
    nbond_prev = 0
    IF(ASSOCIATED(conn_info%bond_a)) nbond_prev = SIZE(conn_info%bond_a)

    IF(iw>0) WRITE(iw,'(T2,A)') 'PSF_INFO| Parsing the NBOND section'
    IF(iw>0) WRITE(iw,'(T2,A,I8)') 'PSF_INFO| Previous number of allocated BOND: ',nbond_prev
    label = '!NBOND'
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF ( .NOT. found ) THEN
       IF(iw>0) WRITE (iw,'(T2,A)')'PSF_INFO| No NBOND section '
       nbond = 0
    ELSE
       CALL parser_get_object  (parser,nbond,error=error)
       IF(iw>0) WRITE(iw,'(T2,A,'//TRIM(c_int)//')')'PSF_INFO| NBOND = ',nbond
       !malloc the memory that we need
       CALL reallocate(conn_info%bond_a,1,nbond_prev+nbond)
       CALL reallocate(conn_info%bond_b,1,nbond_prev+nbond)
       !Read in the bond info
       IF (psf_type==do_conn_psf_u) THEN
          DO ibond=1,nbond,4
             CALL parser_get_next_line(parser,1,error=error)
             index_now=nbond_prev+ibond-1
             READ(parser%input_line,FMT=*,ERR=9)(conn_info%bond_a(index_now+i),&
                                                 conn_info%bond_b(index_now+i),&
                                                 i=1,MIN(4,(nbond-ibond+1)))
          END DO
       ELSE
          DO ibond=1,nbond,4
             CALL parser_get_next_line(parser,1,error=error)
             index_now=nbond_prev+ibond-1
             READ(parser%input_line,FMT='(8'//TRIM(c_int)//')')&
                  (conn_info%bond_a(index_now+i),&
                  conn_info%bond_b(index_now+i),&
                  i=1,MIN(4,(nbond-ibond+1)))
          END DO
       END IF
       IF ( ANY(conn_info%bond_a(nbond_prev+1:)<=0) .OR. &
            ANY(conn_info%bond_a(nbond_prev+1:)>natom) .OR. &
            ANY(conn_info%bond_b(nbond_prev+1:)<=0) .OR. &
            ANY(conn_info%bond_b(nbond_prev+1:)>natom)) THEN
          CPAssert(.FALSE.,cp_failure_level,"topology_read, invalid bond",error,failure)
       END IF
       conn_info%bond_a(nbond_prev+1:)=conn_info%bond_a(nbond_prev+1:)+natom_prev
       conn_info%bond_b(nbond_prev+1:)=conn_info%bond_b(nbond_prev+1:)+natom_prev
    END IF
    !
    ! THETA section
    !
    ntheta_prev = 0
    IF(ASSOCIATED(conn_info%theta_a)) ntheta_prev = SIZE(conn_info%theta_a)

    IF(iw>0) WRITE(iw,'(T2,A)')'PSF_INFO| Parsing the NTHETA section'
    IF(iw>0) WRITE(iw,'(T2,A,I8)') 'PSF_INFO| Previous number of allocated THETA: ',ntheta_prev
    label = '!NTHETA'
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF ( .NOT. found ) THEN
       IF(iw>0) WRITE (iw,'(T2,A)')'PSF_INFO| No NTHETA section '
       ntheta = 0
    ELSE
       CALL parser_get_object  (parser,ntheta,error=error)
       IF(iw>0) WRITE(iw,'(T2,A,'//TRIM(c_int)//')')'PSF_INFO| NTHETA = ',ntheta
       !malloc the memory that we need
       CALL reallocate(conn_info%theta_a,1,ntheta_prev+ntheta)
       CALL reallocate(conn_info%theta_b,1,ntheta_prev+ntheta)
       CALL reallocate(conn_info%theta_c,1,ntheta_prev+ntheta)
       !Read in the bend info
       IF (psf_type==do_conn_psf_u) THEN
          DO itheta=1,ntheta,3
             CALL parser_get_next_line(parser,1,error=error)
             index_now=ntheta_prev+itheta-1
             READ(parser%input_line,FMT=*,ERR=9)(conn_info%theta_a(index_now+i),&
                                                 conn_info%theta_b(index_now+i),&
                                                 conn_info%theta_c(index_now+i),&
                                                 i=1,MIN(3,(ntheta-itheta+1)))
          END DO
       ELSE
          DO itheta=1,ntheta,3
             CALL parser_get_next_line(parser,1,error=error)
             index_now=ntheta_prev+itheta-1
             READ(parser%input_line,FMT='(9'//TRIM(c_int)//')')&
                  (conn_info%theta_a(index_now+i),&
                  conn_info%theta_b(index_now+i),&
                  conn_info%theta_c(index_now+i),&
                  i=1,MIN(3,(ntheta-itheta+1)))
          END DO
       END IF
       conn_info%theta_a(ntheta_prev+1:)=conn_info%theta_a(ntheta_prev+1:)+natom_prev
       conn_info%theta_b(ntheta_prev+1:)=conn_info%theta_b(ntheta_prev+1:)+natom_prev
       conn_info%theta_c(ntheta_prev+1:)=conn_info%theta_c(ntheta_prev+1:)+natom_prev
    END IF
    !
    ! PHI section
    !
    nphi_prev  = 0
    IF(ASSOCIATED(conn_info%phi_a))   nphi_prev = SIZE(conn_info%phi_a)

    IF(iw>0) WRITE(iw,'(T2,A)') 'PSF_INFO| Parsing the NPHI section'
    IF(iw>0) WRITE(iw,'(T2,A,I8)') 'PSF_INFO| Previous number of allocated PHI: ',nphi_prev
    label = '!NPHI'
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF ( .NOT. found ) THEN
       IF(iw>0) WRITE (iw,'(T2,A)') 'PSF_INFO| No NPHI section '
       nphi = 0
    ELSE
       CALL parser_get_object  (parser,nphi,error=error)
       IF(iw>0) WRITE(iw,'(T2,A,'//TRIM(c_int)//')') 'PSF_INFO| NPHI = ',nphi
       !malloc the memory that we need
       CALL reallocate(conn_info%phi_a,1,nphi_prev+nphi)
       CALL reallocate(conn_info%phi_b,1,nphi_prev+nphi)
       CALL reallocate(conn_info%phi_c,1,nphi_prev+nphi)
       CALL reallocate(conn_info%phi_d,1,nphi_prev+nphi)
       !Read in the torsion info
       IF (psf_type==do_conn_psf_u) THEN
          DO iphi=1,nphi,2
             CALL parser_get_next_line(parser,1,error=error)
             index_now=nphi_prev+iphi-1
             READ(parser%input_line,FMT=*,ERR=9)(conn_info%phi_a(index_now+i),&
                                                 conn_info%phi_b(index_now+i),&
                                                 conn_info%phi_c(index_now+i),&
                                                 conn_info%phi_d(index_now+i),&
                                                 i=1,MIN(2,(nphi-iphi+1)))
          END DO
       ELSE
          DO iphi=1,nphi,2
             CALL parser_get_next_line(parser,1,error=error)
             index_now=nphi_prev+iphi-1
             READ(parser%input_line,FMT='(8'//TRIM(c_int)//')')&
                  (conn_info%phi_a(index_now+i),&
                  conn_info%phi_b(index_now+i),&
                  conn_info%phi_c(index_now+i),&
                  conn_info%phi_d(index_now+i),&
                  i=1,MIN(2,(nphi-iphi+1)))
          END DO
       END IF
       conn_info%phi_a(nphi_prev+1:)=conn_info%phi_a(nphi_prev+1:)+natom_prev
       conn_info%phi_b(nphi_prev+1:)=conn_info%phi_b(nphi_prev+1:)+natom_prev
       conn_info%phi_c(nphi_prev+1:)=conn_info%phi_c(nphi_prev+1:)+natom_prev
       conn_info%phi_d(nphi_prev+1:)=conn_info%phi_d(nphi_prev+1:)+natom_prev
    END IF
    !
    ! IMPHI section
    !
    nphi_prev = 0
    IF(ASSOCIATED(conn_info%impr_a)) nphi_prev = SIZE(conn_info%impr_a)

    IF(iw>0) WRITE(iw,'(T2,A)')'PSF_INFO| Parsing the NIMPHI section'
    IF(iw>0) WRITE(iw,'(T2,A,I8)') 'PSF_INFO| Previous number of allocated IMPHI: ',nphi_prev
    label = '!NIMPHI'
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF ( .NOT. found ) THEN
       IF(iw>0) WRITE (iw,'(T2,A)')'PSF_INFO| No NIMPHI section '
       nphi = 0
    ELSE
       CALL parser_get_object  (parser,nphi,error=error)
       IF(iw>0) WRITE(iw,'(T2,A,'//TRIM(c_int)//')')'PSF_INFO| NIMPR = ',nphi
       !malloc the memory that we need
       CALL reallocate(conn_info%impr_a,1,nphi_prev+nphi)
       CALL reallocate(conn_info%impr_b,1,nphi_prev+nphi)
       CALL reallocate(conn_info%impr_c,1,nphi_prev+nphi)
       CALL reallocate(conn_info%impr_d,1,nphi_prev+nphi)
       !Read in the improper torsion info
       IF (psf_type==do_conn_psf_u) THEN
          DO iphi=1,nphi,2
             CALL parser_get_next_line(parser,1,error=error)
             index_now=nphi_prev+iphi-1
             READ(parser%input_line,FMT=*,ERR=9)(conn_info%impr_a(index_now+i),&
                                                 conn_info%impr_b(index_now+i),&
                                                 conn_info%impr_c(index_now+i),&
                                                 conn_info%impr_d(index_now+i),&
                                                 i=1,MIN(2,(nphi-iphi+1)))
          END DO
       ELSE
          DO iphi=1,nphi,2
             CALL parser_get_next_line(parser,1,error=error)
             index_now=nphi_prev+iphi-1
             READ(parser%input_line,FMT='(8'//TRIM(c_int)//')')&
                  (conn_info%impr_a(index_now+i),&
                  conn_info%impr_b(index_now+i),&
                  conn_info%impr_c(index_now+i),&
                  conn_info%impr_d(index_now+i),&
                  i=1,MIN(2,(nphi-iphi+1)))
          END DO
       END IF
       conn_info%impr_a(nphi_prev+1:)=conn_info%impr_a(nphi_prev+1:)+natom_prev
       conn_info%impr_b(nphi_prev+1:)=conn_info%impr_b(nphi_prev+1:)+natom_prev
       conn_info%impr_c(nphi_prev+1:)=conn_info%impr_c(nphi_prev+1:)+natom_prev
       conn_info%impr_d(nphi_prev+1:)=conn_info%impr_d(nphi_prev+1:)+natom_prev
    END IF

    CALL parser_release(parser,error=error)
    CALL timestop(handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/PSF_INFO",error=error)
    RETURN
9   CONTINUE
    ! Print error and exit
    IF (output_unit>0) THEN
       WRITE(output_unit,'(T2,A)')&
            "PSF_INFO| Error while reading PSF using the unformatted PSF reading option!",&
            "PSF_INFO| Try using PSF instead of UPSF."
    END IF

    CALL stop_program(routineN,moduleN,__LINE__,"Error while reading PSF data!")

  END SUBROUTINE read_topology_psf

! *****************************************************************************
!> \brief Post processing of PSF informations
!> \param topology ...
!> \param subsys_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE  psf_post_process(topology,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'psf_post_process', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iatom, ibond, &
                                                ionfo, iw, jatom, N, natom, &
                                                nbond, nonfo, nphi, ntheta, &
                                                stat
    LOGICAL                                  :: failure
    TYPE(array1_list_type), DIMENSION(:), &
      POINTER                                :: ex_bend_list, ex_bond_list
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    failure = .FALSE.
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/PSF_INFO",&
         extension=".subsysLog",error=error)
    CALL timeset(routineN,handle)
    atom_info => topology%atom_info
    conn_info => topology%conn_info
    !
    ! PARA_RES structure
    !
    natom = 0
    nbond = 0
    i     = 0
    IF(ASSOCIATED(atom_info%id_molname)) natom = SIZE(atom_info%id_molname)
    IF (ASSOCIATED(conn_info%bond_a)) nbond = SIZE(conn_info%bond_a)
    IF(ASSOCIATED(conn_info%c_bond_a))    i = SIZE(conn_info%c_bond_a)
    DO ibond=1,nbond
       iatom = conn_info%bond_a(ibond)
       jatom = conn_info%bond_b(ibond)
       IF(topology%para_res) THEN
          IF((atom_info%id_molname(iatom)/=atom_info%id_molname(jatom)).OR.&
               (atom_info%resid(iatom)/=atom_info%resid(jatom)).OR.&
               (atom_info%id_resname(iatom)/=atom_info%id_resname(jatom)))THEN
             IF(iw>0) WRITE(iw,'(T2,A,2I6)') "PSF_INFO| PARA_RES, bond between molecules atom ",&
                  iatom,jatom
             i = i + 1
             CALL reallocate(conn_info%c_bond_a,1,i)
             CALL reallocate(conn_info%c_bond_b,1,i)
             conn_info%c_bond_a(i) = iatom
             conn_info%c_bond_b(i) = jatom
          END IF
       ELSE
          IF(atom_info%id_molname(iatom)/=atom_info%id_molname(jatom)) THEN
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END IF
       END IF
    END DO
    !
    ! UB structure
    !
    ntheta = 0
    IF (ASSOCIATED(conn_info%theta_a)) ntheta = SIZE(conn_info%theta_a)
    CALL reallocate(conn_info%ub_a,1,ntheta)
    CALL reallocate(conn_info%ub_b,1,ntheta)
    CALL reallocate(conn_info%ub_c,1,ntheta)
    conn_info%ub_a(:) = conn_info%theta_a(:)
    conn_info%ub_b(:) = conn_info%theta_b(:)
    conn_info%ub_c(:) = conn_info%theta_c(:)
    !
    ! ONFO structure
    !
    nphi = 0
    nonfo= 0
    IF (ASSOCIATED(conn_info%phi_a)) nphi = SIZE(conn_info%phi_a)
    CALL reallocate(conn_info%onfo_a,1,nphi)
    CALL reallocate(conn_info%onfo_b,1,nphi)
    conn_info%onfo_a(1:) = conn_info%phi_a(1:)
    conn_info%onfo_b(1:) = conn_info%phi_d(1:)
    ! Reorder bonds
    ALLOCATE(ex_bond_list(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO I=1,natom
       ALLOCATE(ex_bond_list(I)%array1(0),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDDO
    N = 0
    IF(ASSOCIATED(conn_info%bond_a)) N = SIZE(conn_info%bond_a)
    CALL reorder_structure(ex_bond_list, conn_info%bond_a, conn_info%bond_b, N, error)
    ! Reorder bends
    ALLOCATE(ex_bend_list(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO I=1,natom
       ALLOCATE(ex_bend_list(I)%array1(0),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDDO
    N = 0
    IF(ASSOCIATED(conn_info%theta_a)) N = SIZE(conn_info%theta_a)
    CALL reorder_structure(ex_bend_list, conn_info%theta_a, conn_info%theta_c, N, error)
    DO ionfo=1, nphi
    ! Check if the torsion is not shared between angles or bonds
       IF ( ANY(ex_bond_list(conn_info%onfo_a(ionfo))%array1==conn_info%onfo_b(ionfo)).OR.&
            ANY(ex_bend_list(conn_info%onfo_a(ionfo))%array1==conn_info%onfo_b(ionfo))) CYCLE
       nonfo = nonfo + 1
       conn_info%onfo_a(nonfo) = conn_info%onfo_a(ionfo)
       conn_info%onfo_b(nonfo) = conn_info%onfo_b(ionfo)
    END DO
    ! deallocate bends
    DO I=1,natom
       DEALLOCATE(ex_bend_list(I)%array1,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDDO
    DEALLOCATE(ex_bend_list,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! deallocate bonds
    DO I=1,natom
       DEALLOCATE(ex_bond_list(I)%array1,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDDO
    DEALLOCATE(ex_bond_list,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! Get unique onfo
    ALLOCATE(ex_bond_list(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO I=1,natom
       ALLOCATE(ex_bond_list(I)%array1(0),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDDO
    N = 0
    IF(ASSOCIATED(conn_info%onfo_a)) N = nonfo
    CALL reorder_structure(ex_bond_list, conn_info%onfo_a, conn_info%onfo_b, N, error)
    nonfo = 0
    DO I=1,natom
       DO ionfo = 1, SIZE(ex_bond_list(I)%array1)
          IF (COUNT(ex_bond_list(I)%array1==ex_bond_list(I)%array1(ionfo))/=1) THEN
             ex_bond_list(I)%array1(ionfo) = 0
          ELSE
             IF (ex_bond_list(I)%array1(ionfo)<=I) CYCLE
             nonfo = nonfo + 1
             conn_info%onfo_a(nonfo) = I
             conn_info%onfo_b(nonfo) = ex_bond_list(I)%array1(ionfo)
          END IF
       END DO
    END DO
    DO I=1,natom
       DEALLOCATE(ex_bond_list(I)%array1,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDDO
    DEALLOCATE(ex_bond_list,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL reallocate(conn_info%onfo_a,1,nonfo)
    CALL reallocate(conn_info%onfo_b,1,nonfo)

    CALL timestop(handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/PSF_INFO",error=error)
  END SUBROUTINE psf_post_process

! *****************************************************************************
!> \brief Input driven modification (IDM) of PSF defined structures
!> \param topology ...
!> \param section ...
!> \param subsys_section ...
!> \param error ...
!> \author Teodoro Laino - Zurich University 04.2007
! *****************************************************************************
  SUBROUTINE  idm_psf(topology,section,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(section_vals_type), POINTER         :: section, subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'idm_psf', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iend, iend1, istart, istart1, item, iw, j, mol_id, &
      n_rep, natom, nbond, nimpr, noe, nphi, ntheta, stat
    INTEGER, DIMENSION(:), POINTER           :: tag_mols, tmp, wrk
    LOGICAL                                  :: explicit, failure
    TYPE(array1_list_type), DIMENSION(:), &
      POINTER                                :: ex_bond_list
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: subsection

    NULLIFY(logger)
    failure = .FALSE.
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/PSF_INFO",&
         extension=".subsysLog",error=error)
    CALL timeset(routineN,handle)
    CALL section_vals_get(section, explicit=explicit, error=error)
    IF (explicit) THEN
       atom_info => topology%atom_info
       conn_info => topology%conn_info
       natom = 0
       IF(ASSOCIATED(atom_info%id_molname)) natom = SIZE(atom_info%id_molname)
       nbond = 0
       IF (ASSOCIATED(conn_info%bond_a)) nbond = SIZE(conn_info%bond_a)
       ntheta = 0
       IF (ASSOCIATED(conn_info%theta_a)) ntheta = SIZE(conn_info%theta_a)
       nphi = 0
       IF (ASSOCIATED(conn_info%phi_a)) nphi = SIZE(conn_info%phi_a)
       nimpr = 0
       IF (ASSOCIATED(conn_info%impr_a)) nimpr = SIZE(conn_info%impr_a)
       ! Any new defined bond
       subsection => section_vals_get_subs_vals(section,"BONDS",error=error)
       CALL section_vals_get(subsection, explicit=explicit, error=error)
       IF (explicit) THEN
          CALL section_vals_val_get(subsection,"_DEFAULT_KEYWORD_",n_rep_val=n_rep, error=error)
          CALL reallocate(conn_info%bond_a,1,n_rep+nbond)
          CALL reallocate(conn_info%bond_b,1,n_rep+nbond)
          DO i = 1, n_rep
             CALL section_vals_val_get(subsection,"_DEFAULT_KEYWORD_",i_rep_val=i, i_vals=tmp, error=error)
             conn_info%bond_a(nbond+i)=tmp(1)
             conn_info%bond_b(nbond+i)=tmp(2)
          END DO
          ! And now modify the molecule name if two molecules have been bridged
          ALLOCATE(ex_bond_list(natom),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(tag_mols(natom),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(wrk(natom),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO j=1,natom
             ALLOCATE(ex_bond_list(j)%array1(0),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDDO
          CALL reorder_structure(ex_bond_list, conn_info%bond_a, conn_info%bond_b, nbond+n_rep, error)
          ! Loop over atoms to possiblyt change molecule name
          tag_mols=-1
          mol_id = 1
          DO i = 1, natom
             IF (tag_mols(i)/=-1) CYCLE
             CALL tag_molecule(tag_mols, ex_bond_list, i, mol_id)
             mol_id = mol_id + 1
          END DO
          mol_id = mol_id - 1
          IF (iw>0) WRITE(iw,'(T2,A,I8)') 'PSF_INFO| Number of molecules detected after merging: ',mol_id
          ! Now simply check about the contiguousness of molecule definition
          CALL sort(tag_mols, natom, wrk)
          item   = tag_mols(1)
          istart = 1
          DO i = 2, natom
             IF (tag_mols(i)==item) CYCLE
             iend = i-1
             noe = iend-istart+1
             istart1 = MINVAL(wrk(istart:iend))
             iend1   = MAXVAL(wrk(istart:iend))
             CPPostcondition(iend1-istart1+1==noe,cp_failure_level,routineP,error,failure)
             atom_info%id_molname(istart1:iend1)=str2id(s2s("MOL"//cp_to_string(item)))
             item = tag_mols(i)
             istart = i
          END DO
          iend = i-1
          noe  = iend-istart+1
          istart1 = MINVAL(wrk(istart:iend))
          iend1   = MAXVAL(wrk(istart:iend))
          CPPostcondition(iend1-istart1+1==noe,cp_failure_level,routineP,error,failure)
          atom_info%id_molname(istart1:iend1)=str2id(s2s("MOL"//cp_to_string(item)))
          ! Deallocate bonds
          DO i=1,natom
             DEALLOCATE(ex_bond_list(i)%array1,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDDO
          DEALLOCATE(ex_bond_list,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(tag_mols,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(wrk,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       ! Any new defined angle
       subsection => section_vals_get_subs_vals(section,"ANGLES",error=error)
       CALL section_vals_get(subsection, explicit=explicit, error=error)
       IF (explicit) THEN
          CALL section_vals_val_get(subsection,"_DEFAULT_KEYWORD_",n_rep_val=n_rep, error=error)
          CALL reallocate(conn_info%theta_a,1,n_rep+ntheta)
          CALL reallocate(conn_info%theta_b,1,n_rep+ntheta)
          CALL reallocate(conn_info%theta_c,1,n_rep+ntheta)
          DO i = 1, n_rep
             CALL section_vals_val_get(subsection,"_DEFAULT_KEYWORD_",i_rep_val=i, i_vals=tmp, error=error)
             conn_info%theta_a(ntheta+i)=tmp(1)
             conn_info%theta_b(ntheta+i)=tmp(2)
             conn_info%theta_c(ntheta+i)=tmp(3)
          END DO
       END IF
       ! Any new defined torsion
       subsection => section_vals_get_subs_vals(section,"TORSIONS",error=error)
       CALL section_vals_get(subsection, explicit=explicit, error=error)
       IF (explicit) THEN
          CALL section_vals_val_get(subsection,"_DEFAULT_KEYWORD_",n_rep_val=n_rep, error=error)
          CALL reallocate(conn_info%phi_a,1,n_rep+nphi)
          CALL reallocate(conn_info%phi_b,1,n_rep+nphi)
          CALL reallocate(conn_info%phi_c,1,n_rep+nphi)
          CALL reallocate(conn_info%phi_d,1,n_rep+nphi)
          DO i = 1, n_rep
             CALL section_vals_val_get(subsection,"_DEFAULT_KEYWORD_",i_rep_val=i, i_vals=tmp, error=error)
             conn_info%phi_a(nphi+i)=tmp(1)
             conn_info%phi_b(nphi+i)=tmp(2)
             conn_info%phi_c(nphi+i)=tmp(3)
             conn_info%phi_d(nphi+i)=tmp(4)
          END DO
       END IF
       ! Any new defined improper
       subsection => section_vals_get_subs_vals(section,"IMPROPERS",error=error)
       CALL section_vals_get(subsection, explicit=explicit, error=error)
       IF (explicit) THEN
          CALL section_vals_val_get(subsection,"_DEFAULT_KEYWORD_",n_rep_val=n_rep, error=error)
          CALL reallocate(conn_info%impr_a,1,n_rep+nimpr)
          CALL reallocate(conn_info%impr_b,1,n_rep+nimpr)
          CALL reallocate(conn_info%impr_c,1,n_rep+nimpr)
          CALL reallocate(conn_info%impr_d,1,n_rep+nimpr)
          DO i = 1, n_rep
             CALL section_vals_val_get(subsection,"_DEFAULT_KEYWORD_",i_rep_val=i, i_vals=tmp, error=error)
             conn_info%impr_a(nimpr+i)=tmp(1)
             conn_info%impr_b(nimpr+i)=tmp(2)
             conn_info%impr_c(nimpr+i)=tmp(3)
             conn_info%impr_d(nimpr+i)=tmp(4)
          END DO
       END IF
    END IF
    CALL timestop(handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/PSF_INFO",error=error)

  END SUBROUTINE idm_psf

! *****************************************************************************
!> \brief Teodoro Laino - 01.2006
!>      Write PSF topology file in the CHARMM31 EXT standard format
!> \param file_unit ...
!> \param topology ...
!> \param subsys_section ...
!> \param force_env_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE write_topology_psf (file_unit,topology,subsys_section,force_env_section,error)
    INTEGER, INTENT(IN)                      :: file_unit
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(section_vals_type), POINTER         :: subsys_section, &
                                                force_env_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_topology_psf', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2*default_string_length)   :: psf_format
    CHARACTER(LEN=default_string_length)     :: c_int, my_tag1, my_tag2, &
                                                my_tag3, record
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: charge_atm
    INTEGER                                  :: handle, i, iw, j, my_index, &
                                                nchg, stat
    LOGICAL                                  :: explicit, failure, ldum
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charge_inp, charges
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key, tmp_section

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    print_key => section_vals_get_subs_vals(subsys_section,"TOPOLOGY%DUMP_PSF",error=error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/PSF_INFO",&
         extension=".subsysLog",error=error)
    CALL timeset(routineN,handle)

    atom_info => topology%atom_info
    conn_info => topology%conn_info

    ! Check for charges.. (need to dump them in the PSF..)
    ALLOCATE(charges(topology%natoms),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    charges = atom_info%atm_charge
    ! Collect charges from Input file..
    NULLIFY(tmp_section)
    tmp_section => section_vals_get_subs_vals(force_env_section,"MM%FORCEFIELD%CHARGE",error=error)
    CALL section_vals_get(tmp_section,explicit=explicit,n_repetition=nchg,error=error)
    IF (explicit) THEN
       ALLOCATE(charge_atm(nchg),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(charge_inp(nchg),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL read_chrg_section(charge_atm,charge_inp,section=tmp_section,start=0,error=error)
       DO j = 1, topology%natoms
          record = id2str(atom_info%id_atmname(j))
          ldum = qmmm_ff_precond_only_qm(record)
          CALL uppercase(record)
          DO i = 1, nchg
             IF (record==charge_atm(i)) THEN
                charges(j) = charge_inp(i)
                EXIT
             END IF
          END DO
       END DO
       DEALLOCATE(charge_atm,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(charge_inp,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    ! fixup for topology output
    DO j = 1, topology%natoms
       IF (charges(j).EQ.-HUGE(0.0_dp)) charges(j)=-99.0_dp
    ENDDO
    record = cp_print_key_generate_filename(logger,print_key,&
         extension=".psf",my_local=.FALSE.,error=error)
    ! build the EXT format
    c_int="I10"
    psf_format = '(I10,T12,A,T21,I0,T30,A,T39,A,T47,A,T53,F10.6,T69,F8.3,T88,I1)'
    IF(iw>0) WRITE(iw,'(T2,A)') &
         "PSF_WRITE| Writing out PSF file with CHARMM31 EXTErnal format: ",TRIM(record)

    WRITE(file_unit,FMT='(A)') "PSF EXT"
    WRITE(file_unit,FMT='(A)') ""
    WRITE(file_unit,FMT='('//TRIM(c_int)//',A)') 1," !NTITLE"
    WRITE(file_unit,FMT='(A)') "   CP2K generated DUMP of connectivity"
    WRITE(file_unit,FMT='(A)') ""

    WRITE(file_unit,FMT='('//TRIM(c_int)//',A)')topology%natoms," !NATOM"
    my_index = 1
    i        = 1
    my_tag1 = id2str(atom_info%id_molname(i))
    my_tag2 = id2str(atom_info%id_resname(i))
    my_tag3 = id2str(atom_info%id_atmname(i))
    ldum = qmmm_ff_precond_only_qm(my_tag1)
    ldum = qmmm_ff_precond_only_qm(my_tag2)
    ldum = qmmm_ff_precond_only_qm(my_tag3)
    WRITE(file_unit,FMT=psf_format)&
         i,&
         TRIM(my_tag1),&
         my_index,&
         TRIM(my_tag2),&
         TRIM(my_tag3),&
         TRIM(my_tag3),&
         charges(i),&
         atom_info%atm_mass(i),&
         0
    DO i=2,topology%natoms
       IF  ((atom_info%map_mol_num(i)/=atom_info%map_mol_num(i-1)).OR.&
            (atom_info%map_mol_res(i)/=atom_info%map_mol_res(i-1))) my_index = my_index + 1
       my_tag1 = id2str(atom_info%id_molname(i))
       my_tag2 = id2str(atom_info%id_resname(i))
       my_tag3 = id2str(atom_info%id_atmname(i))
       ldum = qmmm_ff_precond_only_qm(my_tag1)
       ldum = qmmm_ff_precond_only_qm(my_tag2)
       ldum = qmmm_ff_precond_only_qm(my_tag3)
       WRITE(file_unit,FMT=psf_format)&
            i,&
            TRIM(my_tag1),&
            my_index,&
            TRIM(my_tag2),&
            TRIM(my_tag3),&
            TRIM(my_tag3),&
            charges(i),&
            atom_info%atm_mass(i),&
            0
    END DO
    WRITE(file_unit,FMT='(/)')
    DEALLOCATE(charges,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    WRITE(file_unit,FMT='('//TRIM(c_int)//',A)')SIZE(conn_info%bond_a)," !NBOND"
    DO i=1,SIZE(conn_info%bond_a),4
       j=0
       DO WHILE ((j<4).AND.((i+j)<=SIZE(conn_info%bond_a)))
          WRITE(file_unit,FMT='(2('//TRIM(c_int)//'))',ADVANCE="NO") &
               conn_info%bond_a(i+j),conn_info%bond_b(i+j)
          j=j+1
       END DO
       WRITE(file_unit,FMT='(/)',ADVANCE="NO")
    END DO
    WRITE(file_unit,FMT='(/)')

    WRITE(file_unit,FMT='('//TRIM(c_int)//',A)')SIZE(conn_info%theta_a)," !NTHETA"
    DO i=1,SIZE(conn_info%theta_a),3
       j=0
       DO WHILE ((j<3).AND.((i+j)<=SIZE(conn_info%theta_a)))
          WRITE(file_unit,FMT='(3('//TRIM(c_int)//'))',ADVANCE="NO") &
               conn_info%theta_a(i+j),conn_info%theta_b(i+j),&
               conn_info%theta_c(i+j)
          j=j+1
       END DO
       WRITE(file_unit,FMT='(/)',ADVANCE="NO")
    END DO
    WRITE(file_unit,FMT='(/)')

    WRITE(file_unit,FMT='('//TRIM(c_int)//',A)')SIZE(conn_info%phi_a)," !NPHI"
    DO i=1,SIZE(conn_info%phi_a),2
       j=0
       DO WHILE ((j<2).AND.((i+j)<=SIZE(conn_info%phi_a)))
          WRITE(file_unit,FMT='(4('//TRIM(c_int)//'))',ADVANCE="NO") &
               conn_info%phi_a(i+j),conn_info%phi_b(i+j),&
               conn_info%phi_c(i+j),conn_info%phi_d(i+j)
          j=j+1
       END DO
       WRITE(file_unit,FMT='(/)',ADVANCE="NO")
    END DO
    WRITE(file_unit,FMT='(/)')

    WRITE(file_unit,FMT='('//TRIM(c_int)//',A)')SIZE(conn_info%impr_a)," !NIMPHI"
    DO i=1,SIZE(conn_info%impr_a),2
       j=0
       DO WHILE ((j<2).AND.((i+j)<=SIZE(conn_info%impr_a)))
          WRITE(file_unit,FMT='(4('//TRIM(c_int)//'))',ADVANCE="NO") &
               conn_info%impr_a(i+j),conn_info%impr_b(i+j),&
               conn_info%impr_c(i+j),conn_info%impr_d(i+j)
          j=j+1
       END DO
       WRITE(file_unit,FMT='(/)',ADVANCE="NO")
    END DO
    WRITE(file_unit,FMT='(/)')

    WRITE(file_unit,FMT='('//TRIM(c_int)//',A)') 0," !NDON"
    WRITE(file_unit,FMT='(/)')
    WRITE(file_unit,FMT='('//TRIM(c_int)//',A)') 0," !NACC"
    WRITE(file_unit,FMT='(/)')
    WRITE(file_unit,FMT='('//TRIM(c_int)//',A)') 0," !NNB"
    WRITE(file_unit,FMT='(/)')

    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/PSF_INFO",error=error)
    CALL timestop(handle)

  END SUBROUTINE write_topology_psf

END MODULE topology_psf

