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

  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE parser,                          ONLY: get_next,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label
  USE termination,                     ONLY: stop_memory
  USE topology_types,                  ONLY: topology_parameters_type
  USE topology_util,                   ONLY: sort_bend,&
                                             sort_bond,&
                                             sort_dihe

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
    REAL(dbl), POINTER                       :: atm_charge(:), atm_mass(:)

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
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  !get the temporary memory
  ALLOCATE (map_mol_typ(natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','map_mol_typ',nphi)
  ALLOCATE (map_mol_num(natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','map_mol_num',nphi)
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  !Sort the bonds, bends, and dihes
  CALL sort_bond(nbond,bond_a,bond_b)
  CALL sort_bend(ntheta,theta_a,theta_b,theta_c)
  CALL sort_dihe(nphi,phi_a,phi_b,phi_c,phi_d)
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
    REAL(dbl), POINTER                       :: atm_charge(:), atm_mass(:)

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  iw = globenv%scr

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

END SUBROUTINE write_topology_psf


END MODULE topology_psf

!******************************************************************************
