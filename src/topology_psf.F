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

  USE global_types,                     ONLY : global_environment_type
  USE input_types,                      ONLY : setup_parameters_type
  USE kinds,                            ONLY : dbl
  USE message_passing,                  ONLY : mp_bcast
  USE molecule_types,                   ONLY : init_molecule_type,&
                                               molecule_type
  USE nrutil,                           ONLY : swap
  USE parser,                           ONLY : parser_init, &
                                               parser_end, &
                                               read_line, &
                                               test_next, &
                                               p_error, &
                                               get_next, &
                                               search_label
  USE string_utilities,                 ONLY : uppercase, &
                                               xstring, &
                                               str_search
  USE util,                             ONLY : get_unit
  USE termination,                      ONLY : stop_program, &
                                               stop_memory
  USE topology_util,                    ONLY : azero,&
                                               find_boundary,&
                                               sort_bond,&
                                               sort_bend,&
                                               sort_dihe,&
                                               topology_generate_molecule,&
                                               topology_pack
  USE string_utilities,                 ONLY : integer_to_string


  PRIVATE
  PUBLIC :: read_topology_psf

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

SUBROUTINE read_topology_psf (mol_setup,setup,globenv)

  IMPLICIT NONE

! Arguments
  TYPE ( molecule_type ), DIMENSION ( : ), POINTER :: mol_setup
  TYPE ( setup_parameters_type ), INTENT (IN) :: setup
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Locals
  CHARACTER ( LEN=20 ) :: label,str2
  INTEGER :: iw,ierror,ilen
  INTEGER :: iatom,natom
  CHARACTER (LEN=20), POINTER :: atm_mol_name(:)
  INTEGER, POINTER :: atm_res_num(:)
  CHARACTER (LEN=20), POINTER :: atm_res_name(:)
  CHARACTER (LEN=20), POINTER :: atm_name1(:)
  CHARACTER (LEN=20), POINTER :: atm_name2(:)
  REAL (dbl), POINTER :: atm_charge(:)
  REAL (dbl), POINTER :: atm_mass(:)
  INTEGER, POINTER :: map_mol_typ(:) !map_mol_typ(iatom)= mol_typ number
  INTEGER, POINTER :: map_mol_num(:) !map_mol_num(iatom)= mol number
  INTEGER :: ibond,nbond
  INTEGER, POINTER :: bond_a(:),bond_b(:)
  INTEGER :: itheta,ntheta
  INTEGER, POINTER :: theta_a(:),theta_b(:),theta_c(:)
  INTEGER :: iphi,nphi
  INTEGER, POINTER :: phi_a(:),phi_b(:),phi_c(:),phi_d(:)
  INTEGER :: iunit
  INTEGER :: i,j,k,first,last,search
  INTEGER :: nmol_typ
  INTEGER, POINTER :: checkme(:)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  iw = globenv%scr

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  WRITE(iw,*) 'Parsing the NATOM section'

  label = '!NATOM'
  CALL parser_init(setup%topology_file_name,globenv)
  CALL search_label(label,ierror,ignore_case=.TRUE.)
  IF ( ierror /= 0 ) THEN
    WRITE ( iw, '( A )' ) ' No NATOM section '
    natom = 0
  ELSE
    INQUIRE (FILE=setup%topology_file_name,NUMBER=iunit)
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
  CALL parser_init(setup%topology_file_name,globenv)
  CALL search_label(label,ierror,ignore_case=.TRUE.)
  IF ( ierror /= 0 ) THEN
    WRITE ( iw, '( A )' ) ' No NBOND section '
    nbond = 0
  ELSE
    INQUIRE (FILE=setup%topology_file_name,NUMBER=iunit)
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
  CALL parser_init(setup%topology_file_name,globenv)
  CALL search_label(label,ierror,ignore_case=.TRUE.)
  IF ( ierror /= 0 ) THEN
    WRITE ( iw, '( A )' ) ' No NTHETA section '
    ntheta = 0
  ELSE
    INQUIRE (FILE=setup%topology_file_name,NUMBER=iunit)
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
  CALL parser_init(setup%topology_file_name,globenv)
  CALL search_label(label,ierror,ignore_case=.TRUE.)
  IF ( ierror /= 0 ) THEN
    WRITE ( iw, '( A )' ) ' No NTHETA section '
    nphi = 0
  ELSE
    INQUIRE (FILE=setup%topology_file_name,NUMBER=iunit)
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
  CALL topology_generate_molecule(natom,atm_res_name,nbond,bond_a,bond_b,&
                                  map_mol_typ,map_mol_num)
  CALL topology_pack(bond_a,bond_b,theta_a,theta_b,theta_c,&
                     phi_a,phi_b,phi_c,phi_d,&
                     atm_name2,atm_res_name,atm_mass,atm_charge,&
                     map_mol_typ,map_mol_num,&
                     mol_setup,setup,globenv)

END SUBROUTINE read_topology_psf

END MODULE topology_psf

!******************************************************************************
