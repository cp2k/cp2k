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
  ! XXXXXXXXXXXX Stil needed ?? Joost  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !get the temporary memory
  !ALLOCATE (map_mol_typ(natom),STAT=ierror)
  !  IF(ierror/=0) CALL stop_memory &
  !    ('read_topology_psf','map_mol_typ',nphi)
  !ALLOCATE (map_mol_num(natom),STAT=ierror)
  !  IF(ierror/=0) CALL stop_memory &
  !    ('read_topology_psf','map_mol_num',nphi)
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  !Sort the bonds, bends, and dihes
  !CALL sort_bond(nbond,bond_a,bond_b)
  !CALL sort_bend(ntheta,theta_a,theta_b,theta_c)
  !CALL sort_dihe(nphi,phi_a,phi_b,phi_c,phi_d)
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

SUBROUTINE read_topology_mc (topology,globenv)

! used by Matt to handle the toplogy of MC water....put in this subroutine
! because it closely resembles the psf stuff, except nothing is read from file
! this isn't yet complete, as I need to find out how the atom_info stuff
! corresponds to the names here, taken from read_topology_psf
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
