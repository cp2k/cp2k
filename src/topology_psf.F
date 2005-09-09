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
  USE parser,                          ONLY: get_next,&
                                             parser_end,&
                                             parser_init,&
                                             parser_backspace, &
                                             read_line,&
                                             search_label
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
  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "topology_psf"

  PRIVATE
  PUBLIC :: read_topology_psf,write_topology_psf

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

SUBROUTINE read_topology_psf (topology,globenv)


    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv
    TYPE (atom_info_type),POINTER            :: atom_info
    TYPE (connectivity_info_type),POINTER    :: conn_info

    CHARACTER(LEN=default_string_length)              :: label
    INTEGER                                  :: i, iatom, ibond, ierror, &
                                                ilen, iphi, itheta, iunit, &
                                                iw, natom, nbond, nphi, ntheta
    CHARACTER(len=*), PARAMETER :: routineN = 'read_coordinate_psf', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: handle
    LOGICAL                                  :: PM1, PM2, PM3, PM4
    INTEGER :: jatom

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

  IF(PM1) WRITE(iw,*) "  Entering read_topology_psf"

  CALL write_checkpoint_information("entering "//routineN,globenv)
  CALL timeset(routineN,'I','',handle)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  atom_info => topology%atom_info
  conn_info => topology%conn_info

  IF(PM2) WRITE(iw,*) 'Parsing the NATOM section'

  label = '!NATOM'
  CALL parser_init(topology%conn_file_name,globenv)
  CALL search_label(label,ierror,ignore_case=.TRUE.)
  IF ( ierror /= 0 ) THEN
    IF(PM1) WRITE ( iw, '( A )' ) ' No NATOM section '
    natom = 0
  ELSE
    CALL parser_backspace
    CALL read_line
    CALL get_next (natom)
    IF(PM1) WRITE(iw,*) 'NATOM = ',natom
    !malloc the memory that we need
    IF (ASSOCIATED(atom_info%label_molname)) &
         DEALLOCATE(atom_info%label_molname,STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','dealloc :: atm_mol_name',nphi)
    ALLOCATE(atom_info%label_molname(natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','alloc :: atm_mol_name',nphi)

    IF (ASSOCIATED(atom_info%label_resid)) &
         DEALLOCATE(atom_info%label_resid,STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','dealloc :: atm_res_num',nphi)    
    ALLOCATE(atom_info%label_resid(natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','alloc :: atm_res_num',nphi)

    IF (ASSOCIATED(atom_info%label_resname)) &
         DEALLOCATE(atom_info%label_resname,STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','dealloc :: atm_res_name',nphi)    
    ALLOCATE(atom_info%label_resname(natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','alloc :: atm_res_name',nphi)

    IF (ASSOCIATED(atom_info%label_atmname)) &
         DEALLOCATE(atom_info%label_atmname,STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','dealloc :: atm_name1',nphi)
    ALLOCATE(atom_info%label_atmname(natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','alloc :: atm_name1',nphi)

    IF (ASSOCIATED(atom_info%atm_charge)) &
         DEALLOCATE(atom_info%atm_charge,STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','dealloc :: atm_charge',nphi)    
    ALLOCATE(atom_info%atm_charge(natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','alloc :: atm_charge',nphi)

    IF (ASSOCIATED(atom_info%atm_mass)) &
         DEALLOCATE(atom_info%atm_mass,STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','dealloc :: atm_mass',nphi)    
    ALLOCATE(atom_info%atm_mass(natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','alloc :: atm_mass',nphi)

    !Read in the atom info
    DO iatom=1,natom
      ilen = 20
      CALL read_line
      CALL get_next(i)
      CALL get_next(atom_info%label_molname(iatom),ilen)
      CALL get_next(atom_info%label_resid(iatom))
      CALL get_next(atom_info%label_resname(iatom),ilen)
      CALL get_next(atom_info%label_atmname(iatom),ilen)
      CALL get_next(atom_info%label_atmname(iatom),ilen)
      CALL uppercase ( atom_info%label_atmname(iatom) )
      CALL get_next(atom_info%atm_charge(iatom))
      CALL get_next(atom_info%atm_mass(iatom))
    END DO
  END IF
  CALL parser_end
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  IF(PM1) WRITE(iw,*) 'Parsing the NBOND section'
  label = '!NBOND'
  CALL parser_init(topology%conn_file_name,globenv)
  CALL search_label(label,ierror,ignore_case=.TRUE.)
  IF ( ierror /= 0 ) THEN
    IF(PM1) WRITE ( iw, '( A )' ) ' No NBOND section '
    nbond = 0
  ELSE
    CALL parser_backspace
    CALL read_line
    CALL get_next (nbond)
    IF(PM1) WRITE(iw,*) 'NBOND = ',nbond
    !malloc the memory that we need
    ALLOCATE(conn_info%bond_a(nbond),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','bond_a',nbond)
    ALLOCATE(conn_info%bond_b(nbond),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','bond_b',nbond)
    !Read in the atom info
    DO ibond=1,nbond,4
      CALL read_line
      i=0
      DO WHILE ((i<4).AND.((ibond+i)<=nbond))
        CALL get_next(conn_info%bond_a(ibond+i))
        CALL get_next(conn_info%bond_b(ibond+i))
        i=i+1
      END DO
    END DO
  END IF
  CALL parser_end
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
          IF(PM1) WRITE(iw,*) "      PARA_RES, bond between molecules atom ",&
                              iatom,jatom
          i = i + 1
          CALL reallocate(conn_info%c_bond_a,1,i+1)
          CALL reallocate(conn_info%c_bond_b,1,i+1)
          conn_info%c_bond_a(i) = iatom
          conn_info%c_bond_b(i) = jatom
      END IF
    ELSE
      IF(atom_info%label_molname(iatom)/=atom_info%label_molname(jatom)) THEN
          CALL stop_program ("topology_generate_bond",&
                             "bonds between different molecule types??? ")
      END IF
    END IF
  END DO
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  IF(PM1) WRITE(iw,*) 'Parsing the NTHETA section'
  label = '!NTHETA'
  CALL parser_init(topology%conn_file_name,globenv)
  CALL search_label(label,ierror,ignore_case=.TRUE.)
  IF ( ierror /= 0 ) THEN
    IF(PM1) WRITE ( iw, '( A )' ) ' No NTHETA section '
    ntheta = 0
  ELSE
    CALL parser_backspace
    CALL read_line
    CALL get_next (ntheta)
    IF(PM1) WRITE(iw,*) 'NTHETA = ',ntheta
    !malloc the memory that we need
    ALLOCATE(conn_info%theta_a(ntheta),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','theta_a',ntheta)
    ALLOCATE(conn_info%theta_b(ntheta),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','theta_b',ntheta)
    ALLOCATE(conn_info%theta_c(ntheta),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','theta_c',ntheta)
    !Read in the atom info
    DO itheta=1,ntheta,3
      CALL read_line
      i=0
      DO WHILE ((i<3).AND.((itheta+i)<=ntheta))
        CALL get_next(conn_info%theta_a(itheta+i))
        CALL get_next(conn_info%theta_b(itheta+i))
        CALL get_next(conn_info%theta_c(itheta+i))
        i=i+1
      END DO
    END DO
  END IF
  CALL parser_end
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
    ALLOCATE(conn_info%ub_a(ntheta),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','ub_a',ntheta)
    ALLOCATE(conn_info%ub_b(ntheta),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','ub_b',ntheta)
    ALLOCATE(conn_info%ub_c(ntheta),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','ub_c',ntheta)
    conn_info%ub_a(:) = conn_info%theta_a(:)
    conn_info%ub_b(:) = conn_info%theta_b(:)
    conn_info%ub_c(:) = conn_info%theta_c(:)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  IF(PM1) WRITE(iw,*) 'Parsing the NPHI (ie dihedrals) section'
  label = '!NPHI'
  CALL parser_init(topology%conn_file_name,globenv)
  CALL search_label(label,ierror,ignore_case=.TRUE.)
  IF ( ierror /= 0 ) THEN
    IF(PM1) WRITE ( iw, '( A )' ) ' No NPHI section '
    nphi = 0
  ELSE
    CALL parser_backspace()
    CALL read_line
    CALL get_next (nphi)
    IF(PM1) WRITE(iw,*) 'NPHI = ',nphi
    !malloc the memory that we need
    ALLOCATE(conn_info%phi_a(nphi),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','phi_a',nphi)
    ALLOCATE(conn_info%phi_b(nphi),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','phi_b',nphi)
    ALLOCATE(conn_info%phi_c(nphi),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','phi_c',nphi)
    ALLOCATE(conn_info%phi_d(nphi),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','phi_d',nphi)
    !Read in the atom info
    DO iphi=1,nphi,2
      CALL read_line
      i=0
      DO WHILE ((i<2).AND.((iphi+i)<=nphi))
        CALL get_next(conn_info%phi_a(iphi+i))
        CALL get_next(conn_info%phi_b(iphi+i))
        CALL get_next(conn_info%phi_c(iphi+i))
        CALL get_next(conn_info%phi_d(iphi+i))
        i=i+1
      END DO
    END DO
  END IF
  CALL parser_end
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
    ALLOCATE(conn_info%onfo_a(nphi),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','onfo_a',nphi)
    ALLOCATE(conn_info%onfo_b(nphi),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','onfo_b',nphi)
    conn_info%onfo_a(:) = conn_info%phi_a(:)
    conn_info%onfo_b(:) = conn_info%phi_d(:)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  IF(PM1) WRITE(iw,*) 'Parsing the NIMPHI (ie dihedrals) section'
  label = '!NIMPHI'
  CALL parser_init(topology%conn_file_name,globenv)
  CALL search_label(label,ierror,ignore_case=.TRUE.)
  IF ( ierror /= 0 ) THEN
    IF(PM1) WRITE ( iw, '( A )' ) ' No NIMPHI section '
    nphi = 0
  ELSE
    CALL parser_backspace
    CALL read_line
    CALL get_next (nphi)
    IF(PM1) WRITE(iw,*) 'NIMPR = ',nphi
    !malloc the memory that we need
    ALLOCATE(conn_info%impr_a(nphi),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','impr_a',nphi)
    ALLOCATE(conn_info%impr_b(nphi),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','impr_b',nphi)
    ALLOCATE(conn_info%impr_c(nphi),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','impr_c',nphi)
    ALLOCATE(conn_info%impr_d(nphi),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','impr_d',nphi)
    !Read in the atom info
    DO iphi=1,nphi,2
      CALL read_line
      i=0
      DO WHILE ((i<2).AND.((iphi+i)<=nphi))
        CALL get_next(conn_info%impr_a(iphi+i))
        CALL get_next(conn_info%impr_b(iphi+i))
        CALL get_next(conn_info%impr_c(iphi+i))
        CALL get_next(conn_info%impr_d(iphi+i))
        i=i+1
      END DO
    END DO
  END IF
  CALL parser_end

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  IF(PM1) WRITE(iw,*) "  Exiting read_topology_psf"

  CALL timestop(0.0_dp,handle)
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("exiting "//routineN,globenv)
  END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

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

SUBROUTINE write_topology_psf (topology,globenv)


    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    INTEGER                                  :: i, j, iatom, ibond, ierror, &
                                                ilen, iphi, itheta, iunit, &
                                                iw, natom, nbond, nphi, ntheta

    CHARACTER(len=*), PARAMETER :: routineN = 'write_topology_psf', &
      routineP = moduleN//':'//routineN
    TYPE (atom_info_type),POINTER            :: atom_info
    TYPE (connectivity_info_type),POINTER    :: conn_info
    CHARACTER(LEN=default_string_length)     :: record
    INTEGER                                  :: file_unit,imol,id1,id2,handle
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

  IF(PM1) WRITE(iw,*) "  Entering write_topology_psf"

  CALL write_checkpoint_information("entering "//routineN,globenv)
  CALL timeset(routineN,'I','',handle)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  atom_info => topology%atom_info
  conn_info => topology%conn_info

  record = "topology_dump.psf"

  IF(PM1) WRITE(iw,*) "    Writing out PSF file ",TRIM(record)

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
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  IF(PM1) WRITE(iw,*) "  Exiting write_topology_psf"

  CALL timestop(0.0_dp,handle)
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("exiting "//routineN,globenv)
  END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

END SUBROUTINE write_topology_psf

END MODULE topology_psf

!******************************************************************************
