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
  USE molecule_types,                   ONLY : molecule_type
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
                                               topology_generate_molecule
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
      CALL read_line
      CALL get_next(i)
      CALL get_next(atm_mol_name(iatom),ilen)
      CALL get_next(atm_res_num(iatom))
      CALL get_next(atm_res_name(iatom),ilen)
      CALL get_next(atm_name1(iatom),ilen)
      CALL get_next(atm_name2(iatom),ilen)
      CALL get_next(atm_charge(iatom))
      CALL get_next(atm_mass(iatom))
      !WRITE(iw,*) i,atm_res_name(iatom),atm_res_num(iatom),&
      !            atm_res_name(iatom),atm_name1(iatom),atm_name2(iatom),&
      !            atm_charge(iatom),atm_mass(iatom)
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
        !WRITE(iw,*) "AAAA ",ibond,i,ibond+i
        CALL get_next(bond_a(ibond+i))
        CALL get_next(bond_b(ibond+i))
        !WRITE(iw,*) ibond+i,bond_a(ibond+i),bond_b(ibond+i)
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
        !WRITE(iw,*) itheta+i,theta_a(itheta+i),&
        !            theta_b(itheta+i),theta_c(itheta+i)
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
        !WRITE(iw,*) iphi+i,phi_a(iphi+i),&
        !            phi_b(iphi+i),phi_c(iphi+i),phi_d(iphi+i)
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
  !Sort the bonds
  DO ibond=1,nbond
    IF(bond_a(ibond)>bond_b(ibond)) THEN 
      CALL swap(bond_a(ibond),bond_b(ibond))
    END IF
  END DO
  !Need to update this to quicksort for systems greater then 20000 atoms
  DO i=1,nbond
    DO j=1,nbond
      IF(bond_a(i)<bond_a(j)) THEN
        CALL swap(bond_a(i),bond_a(j))
        CALL swap(bond_b(i),bond_b(j))
      END IF
    END DO
  END DO
  i=1
  DO WHILE (i<nbond)
    search = bond_a(i)
    CALL find_boundary(bond_a,nbond,first,last,search)
    !WRITE(iw,*) "find boundary ",search,i,first,last
    IF(.NOT.(first==last)) THEN
      DO j=first,last
        DO k=first,last
          IF(bond_b(j)<bond_b(k)) THEN
            CALL swap(bond_a(j),bond_a(k))
            CALL swap(bond_b(j),bond_b(k))
          END IF
        END DO
      END DO
      i=last+1
    ELSE
      i=i+1
    END IF
  END DO
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  !Sort the angles
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  !Sort the dihedrals
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!Everything down here should be generic----------------------------------------
  WRITE(iw,*) 'calling topology_generate_molecule'
  CALL topology_generate_molecule(natom,atm_res_name,nbond,bond_a,bond_b,&
                                  map_mol_typ,map_mol_num)
!------------------------------------------------------------------------------
  nmol_typ = map_mol_typ(natom)
  IF(nmol_typ < 0) THEN
    CALL stop_program ( "read_topology_psf", "no moleculular type defined")
  ELSE
    !Allocate space for the molecules
    ALLOCATE(mol_setup(nmol_typ),STAT=ierror)
    IF(ierror/=0) CALL stop_memory("read_topology_psf","mol_setup",nmol_typ)
  END IF
!------------------------------------------------------------------------------
  DO i=1,nmol_typ
    CALL find_boundary(map_mol_typ,natom,first,last,i)
    mol_setup(i)%molpar%name = atm_res_name(first)
    mol_setup(i)%num_mol=MAXVAL(map_mol_num)
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !set nhcopt
    mol_setup(i)%nhcopt = "GLOBAL"
    !--------------------------------------------------------------------------
    !Figure out the number of atoms and atom types-----------------------------
    j = last - first + 1
    j = j / map_mol_num(last)
    mol_setup(i)%molpar%natom=j

    mol_setup(i)%molpar%natom_type=1
    DO j=2,mol_setup(i)%molpar%natom
      IF(.NOT.(atm_name2(j-1)==atm_name2(j))) THEN
        mol_setup(i)%molpar%natom_type = mol_setup(i)%molpar%natom_type + 1
      END IF
    END DO
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !get the memory first
    ALLOCATE(mol_setup(i)%molpar%aname(mol_setup(i)%molpar%natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','aname',mol_setup(i)%molpar%natom)
    ALLOCATE(mol_setup(i)%molpar%alabel(mol_setup(i)%molpar%natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','alabel',mol_setup(i)%molpar%natom)
    ALLOCATE(mol_setup(i)%molpar%aweight(mol_setup(i)%molpar%natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','aweight',mol_setup(i)%molpar%natom)
    ALLOCATE(mol_setup(i)%molpar%acharge(mol_setup(i)%molpar%natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','acharge',mol_setup(i)%molpar%natom)
    DO j=1,mol_setup(i)%molpar%natom
      !Stuff the aname,alabel,aweight,acharge
      mol_setup(i)%molpar%aname(j)   = atm_name2(first+j-1)
      mol_setup(i)%molpar%aweight(j) = atm_mass(first+j-1)
      mol_setup(i)%molpar%acharge(j) = atm_charge(first+j-1)
    END DO
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Stuff the bond info
    ALLOCATE(checkme(nbond),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','checkme',nbond)
    CALL azero(checkme,nbond)
    DO j=1,natom
      IF ((map_mol_num(j)==1).AND.(map_mol_typ(j)==i)) THEN
        DO k=1,nbond
          IF ((bond_a(k)==j).OR.(bond_b(k)==j)) THEN
            checkme(k)=1
          END IF
        END DO
      END IF
    END DO
    mol_setup(i)%molpar%nbonds = SUM(checkme)
    WRITE(iw,*) "moltype = ",i,"nbonds = ",mol_setup(i)%molpar%nbonds
    ALLOCATE(mol_setup(i)%molpar%bonds(3,mol_setup(i)%molpar%nbonds),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','mol_setup(i)%molpar%bonds',nbond)
    k=0
    DO j=1,nbond
      IF (checkme(j)==1) THEN
        k=k+1
        !-------------------------------------------------------------------
        !Check and make sure
        mol_setup(i)%molpar%bonds(1,k) = bond_a(j)
        mol_setup(i)%molpar%bonds(2,k) = bond_b(j)
        mol_setup(i)%molpar%bonds(3,k) = 0
        !-------------------------------------------------------------------
      END IF
    END DO
    IF (k /= mol_setup(i)%molpar%nbonds) THEN
      !ERROR STOP PROGRAM
      CALL stop_program ( "read_topology_psf", "nbonds counting error")
    END IF
    DEALLOCATE(checkme,STAT=ierror)
    IF (ierror /= 0) CALL stop_memory ('read_molecule_section','checkme')
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Stuff the theata info
    ALLOCATE(checkme(ntheta),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','checkme',ntheta)
    CALL azero(checkme,ntheta)
    DO j=1,natom
      IF ((map_mol_num(j)==1).AND.(map_mol_typ(j)==i)) THEN
        DO k=1,ntheta
          IF ((theta_a(k)==j).OR.(theta_b(k)==j).OR.(theta_c(k)==j)) THEN
            checkme(k)=1
          END IF
        END DO
      END IF
    END DO
    mol_setup(i)%molpar%nbends = SUM(checkme)
    WRITE(iw,*) "moltype = ",i,"nbends = ",mol_setup(i)%molpar%nbends
    ALLOCATE(mol_setup(i)%molpar%bends(4,mol_setup(i)%molpar%nbends))
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','mol_setup(i)%molpar%bends',&
       mol_setup(i)%molpar%nbends)
    k=0
    DO j=1,ntheta
      IF (checkme(j)==1) THEN
        k=k+1
        !-------------------------------------------------------------------
        !Check and make sure
        mol_setup(i)%molpar%bends(1,k) = theta_a(j)
        mol_setup(i)%molpar%bends(2,k) = theta_b(j)
        mol_setup(i)%molpar%bends(3,k) = theta_c(j)
        mol_setup(i)%molpar%bends(4,k) = 0
        !-------------------------------------------------------------------
      END IF
    END DO
    IF (k /= mol_setup(i)%molpar%nbends) THEN
      !ERROR STOP PROGRAM
      CALL stop_program ( "read_topology_psf", "nbends counting error")
    END IF
    DEALLOCATE(checkme,STAT=ierror)
    IF (ierror /= 0) CALL stop_memory ('read_molecule_section','checkme')
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Stuff the torsion info
    ALLOCATE(checkme(nphi),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('read_topology_psf','checkme',nphi)
    CALL azero(checkme,nphi)
    DO j=1,natom
      IF ((map_mol_num(j)==1).AND.(map_mol_typ(j)==i)) THEN
        DO k=1,nphi
          IF ((phi_a(k)==j).OR.(phi_b(k)==j).OR.(phi_c(k)==j)&
              .OR.(phi_d(k)==j)) THEN
            checkme(k)=1
          END IF
        END DO
      END IF
    END DO
    mol_setup(i)%molpar%ntorsions = SUM(checkme)
    WRITE(iw,*) "moltype = ",i,"ntorsions = ",mol_setup(i)%molpar%ntorsions
    ALLOCATE(mol_setup(i)%molpar%torsions(5,mol_setup(i)%molpar%ntorsions))
    k=0
    DO j=1,nphi
      IF (checkme(j)==1) THEN
        k=k+1
        !-------------------------------------------------------------------
        !Check and make sure
        mol_setup(i)%molpar%torsions(1,k) = phi_a(j)
        mol_setup(i)%molpar%torsions(2,k) = phi_b(j)
        mol_setup(i)%molpar%torsions(3,k) = phi_c(j)
        mol_setup(i)%molpar%torsions(4,k) = phi_d(j)
        mol_setup(i)%molpar%torsions(5,k) = 0
        !-------------------------------------------------------------------
      END IF
    END DO
    IF (k /= mol_setup(i)%molpar%ntorsions) THEN
      !ERROR STOP PROGRAM
      CALL stop_program ( "read_topology_psf", "ntorsions counting error")
    END IF
    DEALLOCATE(checkme,STAT=ierror)
    IF (ierror /= 0) CALL stop_memory ('read_molecule_section','checkme')
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Stuff the ndc info
    mol_setup(i)%molpar%ndcon = 0
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Stuff the n3x3c info
    mol_setup(i)%molpar%n3x3con = 0
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
  END DO

  DEALLOCATE(atm_mol_name,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','atm_mol_name')
  DEALLOCATE(atm_res_num,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','atm_res_num')
  DEALLOCATE(atm_res_name,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','atm_res_name')
  DEALLOCATE(atm_name1,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','atm_name1')
  DEALLOCATE(atm_name2,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','atm_name2')
  DEALLOCATE(atm_charge,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','atm_charge')
  DEALLOCATE(atm_mass,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','atm_mass')
  DEALLOCATE(map_mol_typ,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','map_mol_typ')
  DEALLOCATE(map_mol_num,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','map_mol_num')
  DEALLOCATE(bond_a,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','bond_a')
  DEALLOCATE(bond_b,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','bond_b')
  DEALLOCATE(theta_a,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','theta_a')
  DEALLOCATE(theta_b,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','theta_b')
  DEALLOCATE(theta_c,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','theta_c')
  DEALLOCATE(phi_a,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','phi_a')
  DEALLOCATE(phi_b,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','phi_b')
  DEALLOCATE(phi_c,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','phi_c')
  DEALLOCATE(phi_d,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','phi_d')

  !CALL stop_program ( "read_topology_psf", "checkme")

END SUBROUTINE read_topology_psf

END MODULE topology_psf

!******************************************************************************
