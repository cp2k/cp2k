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
  USE parser,                          ONLY: get_next,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label
  USE qs_parser,                       ONLY: close_file,&
                                             open_file
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: topology_parameters_type,&
                                             atom_info_type,&
                                             connectivity_info_type

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "topology_psf"

  PRIVATE
  PUBLIC :: read_topology_psf,write_topology_psf,read_topology_mc

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
!!   AUTHOR
!!     IKUO 08.01.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE read_topology_psf (topology,globenv)


    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=20)                        :: label
    CHARACTER(LEN=20), POINTER               :: atm_mol_name(:), &
                                                atm_name1(:), atm_name2(:), &
                                                atm_res_name(:)
    INTEGER                                  :: i, iatom, ibond, ierror, &
                                                ilen, iphi, itheta, iunit, &
                                                iw, natom, nbond, nphi, ntheta
    INTEGER, POINTER :: atm_res_num(:), bond_a(:), bond_b(:), map_mol_num(:), &
      map_mol_typ(:), phi_a(:), phi_b(:), phi_c(:), phi_d(:), theta_a(:), &
      theta_b(:), theta_c(:)
    REAL(KIND=dp), POINTER                       :: atm_charge(:), atm_mass(:)

!map_mol_typ(iatom)= mol_typ number
!map_mol_num(iatom)= mol number
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  iw = globenv%scr

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  WRITE(iw,*) 'Parsing the NATOM section'

  label = '!NATOM'
  CALL parser_init(topology%conn_file_name,globenv)
  CALL search_label(label,ierror,ignore_case=.TRUE.)
  IF ( ierror /= 0 ) THEN
    WRITE ( iw, '( A )' ) ' No NATOM section '
    natom = 0
  ELSE
    INQUIRE (FILE=topology%conn_file_name,NUMBER=iunit)
    BACKSPACE iunit 
    CALL read_line
    CALL get_next (natom)
    WRITE(iw,*) 'NATOM = ',natom
    !malloc the memory that we need
    ALLOCATE(atm_mol_name(natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','atm_mol_name',nphi)
    ALLOCATE(atm_res_num(natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','atm_res_num',nphi)
    ALLOCATE(atm_res_name(natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','atm_res_name',nphi)
    ALLOCATE(atm_name1(natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','atm_name1',nphi)
    ALLOCATE(atm_name2(natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','atm_name2',nphi)
    ALLOCATE(atm_charge(natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','atm_charge',nphi)
    ALLOCATE(atm_mass(natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','atm_mass',nphi)
    !Read in the atom info
    DO iatom=1,natom
      ilen = 20
      CALL read_line
      CALL get_next(i)
      CALL get_next(atm_mol_name(iatom),ilen)
      CALL get_next(atm_res_num(iatom))
      CALL get_next(atm_res_name(iatom),ilen)
      CALL get_next(atm_name1(iatom),ilen)
      CALL get_next(atm_name2(iatom),ilen)
      CALL get_next(atm_charge(iatom))
      CALL get_next(atm_mass(iatom))
    END DO
  END IF
  CALL parser_end
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  WRITE(iw,*) 'Parsing the NBOND section'
  label = '!NBOND: bonds'
  CALL parser_init(topology%conn_file_name,globenv)
  CALL search_label(label,ierror,ignore_case=.TRUE.)
  IF ( ierror /= 0 ) THEN
    WRITE ( iw, '( A )' ) ' No NBOND section '
    nbond = 0
  ELSE
    INQUIRE (FILE=topology%conn_file_name,NUMBER=iunit)
    BACKSPACE iunit 
    CALL read_line
    CALL get_next (nbond)
    WRITE(iw,*) 'NBOND = ',nbond
    !malloc the memory that we need
    ALLOCATE(bond_a(nbond),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','bond_a',nphi)
    ALLOCATE(bond_b(nbond),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','bond_b',nphi)
    !Read in the atom info
    DO ibond=1,nbond,4
      CALL read_line
      i=0
      DO WHILE ((i<4).AND.((ibond+i)<=nbond))
        CALL get_next(bond_a(ibond+i))
        CALL get_next(bond_b(ibond+i))
        i=i+1
      END DO
    END DO
  END IF
  CALL parser_end
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  WRITE(iw,*) 'Parsing the NTHETA section'
  label = '!NTHETA: angles'
  CALL parser_init(topology%conn_file_name,globenv)
  CALL search_label(label,ierror,ignore_case=.TRUE.)
  IF ( ierror /= 0 ) THEN
    WRITE ( iw, '( A )' ) ' No NTHETA section '
    ntheta = 0
  ELSE
    INQUIRE (FILE=topology%conn_file_name,NUMBER=iunit)
    BACKSPACE iunit 
    CALL read_line
    CALL get_next (ntheta)
    WRITE(iw,*) 'NTHETA = ',ntheta
    !malloc the memory that we need
    ALLOCATE(theta_a(ntheta),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','theta_a',nphi)
    ALLOCATE(theta_b(ntheta),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','theta_b',nphi)
    ALLOCATE(theta_c(ntheta),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','theta_c',nphi)
    !Read in the atom info
    DO itheta=1,ntheta,3
      CALL read_line
      i=0
      DO WHILE ((i<3).AND.((itheta+i)<=ntheta))
        CALL get_next(theta_a(itheta+i))
        CALL get_next(theta_b(itheta+i))
        CALL get_next(theta_c(itheta+i))
        i=i+1
      END DO
    END DO
  END IF
  CALL parser_end
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  WRITE(iw,*) 'Parsing the NPHI (ie dihedrals) section'
  label = '!NPHI: dihedrals'
  CALL parser_init(topology%conn_file_name,globenv)
  CALL search_label(label,ierror,ignore_case=.TRUE.)
  IF ( ierror /= 0 ) THEN
    WRITE ( iw, '( A )' ) ' No NTHETA section '
    nphi = 0
  ELSE
    INQUIRE (FILE=topology%conn_file_name,NUMBER=iunit)
    BACKSPACE iunit 
    CALL read_line
    CALL get_next (nphi)
    WRITE(iw,*) 'NPHI = ',nphi
    !malloc the memory that we need
    ALLOCATE(phi_a(nphi),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','phi_a',nphi)
    ALLOCATE(phi_b(nphi),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','phi_b',nphi)
    ALLOCATE(phi_c(nphi),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','phi_c',nphi)
    ALLOCATE(phi_d(nphi),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','phi_d',nphi)
    !Read in the atom info
    DO iphi=1,nphi,2
      CALL read_line
      i=0
      DO WHILE ((i<2).AND.((iphi+i)<=nphi))
        CALL get_next(phi_a(iphi+i))
        CALL get_next(phi_b(iphi+i))
        CALL get_next(phi_c(iphi+i))
        CALL get_next(phi_d(iphi+i))
        i=i+1
      END DO
    END DO
  END IF
  CALL parser_end
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!WHEN FIST IS FULLY FUNCTIONAL, JUST READ THIS STUFF IN------------------------
! WRITE(iw,*) 'Parsing the NIMPHI (ie impropers) section'
! WRITE(iw,*) 'Parsing the NDON (ie donors) section'
! WRITE(iw,*) 'Parsing the NACC (ie acceptors) section'
! WRITE(iw,*) 'Parsing the NNB section'
! WRITE(iw,*) 'Parsing the NGRP section'
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!Everything down here should be generic----------------------------------------

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
    INTEGER, POINTER :: atm_res_num(:), bond_a(:), bond_b(:), map_mol_num(:), &
      map_mol_typ(:), phi_a(:), phi_b(:), phi_c(:), phi_d(:), theta_a(:), &
      theta_b(:), theta_c(:)
    REAL(KIND=dp), POINTER                       :: atm_charge(:), atm_mass(:)

    CHARACTER(len=*), PARAMETER :: routineN = 'read_coordinate_psf', &
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

  IF(PM1) WRITE(iw,*) "  Entering write_coordinate_pdb"

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
    WRITE(file_unit,'(I7,1X,A5,I7,1X,A5,A5,A5,F10.5,F10.5,I7)') &
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
      WRITE(file_unit,'(I7,I7)',ADVANCE="NO") &
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
      WRITE(file_unit,'(I7,I7,I7)',ADVANCE="NO") &
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
      WRITE(file_unit,'(I7,I7,I7,I7)',ADVANCE="NO") &
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
      WRITE(file_unit,'(I7,I7,I7,I7)',ADVANCE="NO") &
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
  
  
  !line = ''
  !line(1:6) = "TITLE  "
  !line(11:70) = "CP2K generated DUMP of coordinates"
  !WRITE(file_unit,*) TRIM(line)


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  IF(PM1) WRITE(iw,*) "  Exiting write_coordinate_pdb"

  CALL timestop(0.0_dp,handle)
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("exiting "//routineN,globenv)
  END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

END SUBROUTINE write_topology_psf

SUBROUTINE read_topology_mc (topology,globenv)

! used by Matt to handle the toplogy of MC water....put in this subroutine
! because it closely resembles the psf stuff, except nothing is read from file
! this isn't yet complete, as I need to find out how the atom_info stuff
! corresponds to the names here, taken from read_topology_psf
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=20), POINTER               :: atm_mol_name(:), &
                                                atm_name1(:), atm_name2(:), &
                                                atm_res_name(:)
    INTEGER                                  :: i, iatom, ibond, ierror, &
                                                ilen, iphi, itheta, iunit, &
                                                iw, natom, nbond, nphi, ntheta, istat
    INTEGER, POINTER :: atm_res_num(:), bond_a(:), bond_b(:), map_mol_num(:), &
      map_mol_typ(:), phi_a(:), phi_b(:), phi_c(:), phi_d(:), theta_a(:), &
      theta_b(:), theta_c(:)
    REAL(KIND=dp), POINTER                       :: atm_charge(:), atm_mass(:)
    LOGICAL                                  :: PL0,PL1,PL2,PL3,PL4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  iw = globenv%scr
  PL0=.FALSE.
  PL1=.FALSE.
  PL2=.FALSE.
  PL3=.FALSE.
  PL4=.FALSE.
  IF ( globenv%print_level .GE. 0 .AND. globenv%ionode) PL0 = .TRUE.
  IF ( globenv%print_level .GE. 1 .AND. globenv%ionode) PL1 = .TRUE.
  IF ( globenv%print_level .GE. 2 .AND. globenv%ionode) PL2 = .TRUE.
  IF ( globenv%print_level .GE. 3 .AND. globenv%ionode) PL3 = .TRUE.
  IF ( globenv%print_level .GE. 4 .AND. globenv%ionode) PL4 = .TRUE.
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  IF(PL1) WRITE(iw,*) '    Generating the topology in the MC water way'

  natom=topology%natoms
  IF(PL1) WRITE(iw,*) '    NATOM = ',natom
  !malloc the memory that we need
  ALLOCATE(atm_mol_name(natom),STAT=ierror)
  IF(ierror/=0) CALL stop_memory &
    ('read_topology_mc','atm_mol_name',natom)
  ALLOCATE(atm_res_num(natom),STAT=ierror)
  IF(ierror/=0) CALL stop_memory &
    ('read_topology_mc','atm_res_num',natom)
  ALLOCATE(atm_res_name(natom),STAT=ierror)
  IF(ierror/=0) CALL stop_memory &
    ('read_topology_mc','atm_res_name',natom)
  ALLOCATE(atm_name1(natom),STAT=ierror)
  IF(ierror/=0) CALL stop_memory &
    ('read_topology_mc','atm_name1',natom)
  ALLOCATE(atm_name2(natom),STAT=ierror)
  IF(ierror/=0) CALL stop_memory &
    ('read_topology_mc','atm_name2',natom)
  ALLOCATE(atm_charge(natom),STAT=ierror)
  IF(ierror/=0) CALL stop_memory &
    ('read_topology_mc','atm_charge',natom)
  ALLOCATE(atm_mass(natom),STAT=ierror)
  IF(ierror/=0) CALL stop_memory &
    ('read_topology_mc','atm_mass',natom)

  ! generate the atom info
  DO iatom=1,natom,3

! first oxygen
    atm_mol_name(iatom) = 'Water'
    atm_res_num(iatom) = (iatom+2)/3
    atm_res_name(iatom) = 'QSMC'
    atm_name1(iatom) = 'O'
    atm_name2(iatom) = 'O'
    atm_charge(iatom) = 0.0d0
    atm_mass(iatom) = 15.9994d0

! now the two hydrogens
    atm_mol_name(iatom+1) = 'Water'
    atm_res_num(iatom+1) = (iatom+2)/3
      atm_res_name(iatom+1) = 'QSMC'
      atm_name1(iatom+1) = 'H'
      atm_name2(iatom+1) = 'H'
      atm_charge(iatom+1) = 0.0d0
      atm_mass(iatom+1) = 1.008d0

      atm_mol_name(iatom+2) = 'Water'
      atm_res_num(iatom+2) = (iatom+2)/3
      atm_res_name(iatom+2) = 'QSMC'
      atm_name1(iatom+2) = 'H'
      atm_name2(iatom+2) = 'H'
      atm_charge(iatom+2) = 0.0d0
      atm_mass(iatom+2) = 1.008d0

  END DO

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  nbond=natom/3*2
  IF(PL1) WRITE(iw,*) '    NBOND = ',nbond
  !malloc the memory that we need
  ALLOCATE(bond_a(nbond),STAT=ierror)
  IF(ierror/=0) CALL stop_memory &
    ('read_topology_mc','bond_a',nphi)
  ALLOCATE(bond_b(nbond),STAT=ierror)
  IF(ierror/=0) CALL stop_memory &
    ('read_topology_mc','bond_b',nphi)
  ! make the bond info
  ibond=1
  DO iatom=1,natom,3
      bond_a(ibond)=iatom
      bond_b(ibond)=iatom+1
      bond_a(ibond+1)=iatom
      bond_b(ibond+1)=iatom+2

      ibond=ibond+2
  END DO

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  ntheta=natom/3
  IF(PL1) WRITE(iw,*) '    NTHETA = ',ntheta
  !malloc the memory that we need
  ALLOCATE(theta_a(ntheta),STAT=ierror)
  IF(ierror/=0) CALL stop_memory &
    ('read_topology_mc','theta_a',ntheta)
  ALLOCATE(theta_b(ntheta),STAT=ierror)
  IF(ierror/=0) CALL stop_memory &
    ('read_topology_mc','theta_b',ntheta)
  ALLOCATE(theta_c(ntheta),STAT=ierror)
  IF(ierror/=0) CALL stop_memory &
    ('read_topology_mc','theta_c',ntheta)
  !Read in the atom info
  itheta=1
  DO iatom=1,natom,3
      theta_a(itheta)=iatom+1
      theta_b(itheta)=iatom
      theta_c(itheta)=iatom+2
      itheta=itheta+1
  END DO

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! no dihedrals in water
  nphi = 0
  IF(PL1) WRITE(iw,*) '    NPHI = ',nphi

  !Sort the bonds, bends, and dihes
  !CALL sort_bond(nbond,bond_a,bond_b)
  !CALL sort_bend(ntheta,theta_a,theta_b,theta_c)
  !CALL sort_dihe(nphi,phi_a,phi_b,phi_c,phi_d)

  !Copy the stuff from temporary array into the more permanent one
  ALLOCATE(topology%conn_info%bond_a(nbond),STAT=istat)
  IF(istat/=0) CALL stop_memory('read_topology_mc','bond_a',nbond)
  ALLOCATE(topology%conn_info%bond_b(nbond),STAT=istat)
  IF(istat/=0) CALL stop_memory('read_topology_mc','bond_b',nbond)
  DO ibond=1,nbond
    topology%conn_info%bond_a(ibond) = bond_a(ibond)
    topology%conn_info%bond_b(ibond) = bond_b(ibond)
  END DO
  ALLOCATE(topology%conn_info%theta_a(ntheta),STAT=istat)
  IF(istat/=0) CALL stop_memory('read_topology_mc','theta_a',ntheta)
  ALLOCATE(topology%conn_info%theta_b(ntheta),STAT=istat)
  IF(istat/=0) CALL stop_memory('read_topology_mc','theta_b',ntheta)
  ALLOCATE(topology%conn_info%theta_c(ntheta),STAT=istat)
  IF(istat/=0) CALL stop_memory('read_topology_mc','theta_c',ntheta)
  DO itheta=1,ntheta
    topology%conn_info%theta_a(itheta) = theta_a(itheta)
    topology%conn_info%theta_b(itheta) = theta_b(itheta)
    topology%conn_info%theta_c(itheta) = theta_c(itheta)
  END DO
! and torsional stuff, even though we have none...just so it's associated
  ALLOCATE(topology%conn_info%phi_a(0),STAT=istat)
  IF(istat/=0) CALL stop_memory('read_topology_mc','phi_a',0)
  ALLOCATE(topology%conn_info%phi_b(0),STAT=istat)
  IF(istat/=0) CALL stop_memory('read_topology_mc','phi_b',0)
  ALLOCATE(topology%conn_info%phi_c(0),STAT=istat)
  IF(istat/=0) CALL stop_memory('read_topology_mc','phi_c',0)
  ALLOCATE(topology%conn_info%phi_d(0),STAT=istat)
  IF(istat/=0) CALL stop_memory('read_topology_mc','phi_d',0)

!------------------------------------------------------------------------------

! now take care of the atom_info
  ALLOCATE (topology%atom_info%map_mol_typ(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_topology_mc','map_mol_typ')
  ALLOCATE (topology%atom_info%map_mol_num(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_topology_mc','map_mol_num')

! we only have one molecule type
  topology%atom_info%map_mol_typ(:) = 1

! the atoms in a molecule are all next to each other
  topology%atom_info%map_mol_num(:) = 0
  DO iatom=1,natom,3
      topology%atom_info%map_mol_num(iatom)=(iatom+2)/3
      topology%atom_info%map_mol_num(iatom+1)=(iatom+2)/3
      topology%atom_info%map_mol_num(iatom+2)=(iatom+2)/3
  ENDDO


! deallocate all the temporary stuff
  DEALLOCATE(bond_a,STAT=istat)
  IF(istat/=0) CALL stop_memory('read_topology_mc','bond_a')
  DEALLOCATE(bond_b,STAT=istat)
  IF(istat/=0) CALL stop_memory('read_topology_mc','bond_b')
  DEALLOCATE(theta_a,STAT=istat)
  IF(istat/=0) CALL stop_memory('read_topology_mc','theta_a')
  DEALLOCATE(theta_b,STAT=istat)
  IF(istat/=0) CALL stop_memory('read_topology_mc','theta_b')
  DEALLOCATE(theta_c,STAT=istat)
  IF(istat/=0) CALL stop_memory('read_topology_mc','theta_c')
  !dellocate some more stuff...this never gets used outside this routine
  DEALLOCATE(atm_mol_name,STAT=ierror)
  IF(ierror/=0) CALL stop_memory('read_topology_mc','atm_mol_name')
  DEALLOCATE(atm_res_num,STAT=ierror)
  IF(ierror/=0) CALL stop_memory('read_topology_mc','atm_res_num')
  DEALLOCATE(atm_res_name,STAT=ierror)
  IF(ierror/=0) CALL stop_memory('read_topology_mc','atm_res_name')
  DEALLOCATE(atm_name1,STAT=ierror)
  IF(ierror/=0) CALL stop_memory('read_topology_mc','atm_name1')
  DEALLOCATE(atm_name2,STAT=ierror)
  IF(ierror/=0) CALL stop_memory('read_topology_mc','atm_name2')
  DEALLOCATE(atm_charge,STAT=ierror)
  IF(ierror/=0) CALL stop_memory('read_topology_mc','atm_charge')
  DEALLOCATE(atm_mass,STAT=ierror)
  IF(ierror/=0) CALL stop_memory('read_topology_mc','atm_mass')
END SUBROUTINE read_topology_mc

END MODULE topology_psf

!******************************************************************************
