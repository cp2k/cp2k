!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/topology_util [1.0] *
!!
!!   NAME
!!     topology_util
!!
!!   FUNCTION
!!     Collection of subroutine needed for topology related things
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!     
!!
!!   SOURCE
!******************************************************************************

MODULE topology_util

  USE atomic_kind_types,               ONLY: allocate_atomic_kind_set,&
                                             atomic_kind_type,&
                                             get_atomic_kind,&
                                             init_atomic_kind_set,&
                                             set_atomic_kind
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE external_potential_types,        ONLY: elp_potential_type,&
                                             set_potential
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl,&
                                             default_string_length
  USE memory_utilities,                ONLY: reallocate
  USE molecule_kind_types,             ONLY: allocate_molecule_kind_set,&
                                             atom_type,&
                                             bend_type,&
                                             bond_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type,&
                                             set_molecule_kind,&
                                             torsion_type
  USE molecule_types_new,              ONLY: allocate_molecule_set,&
                                             local_molecule_type,&
                                             molecule_type,&
                                             set_molecule,&
                                             set_molecule_set
  USE nrutil,                          ONLY: swap
  USE particle_types,                  ONLY: allocate_particle_set,&
                                             particle_type
  USE periodic_table,                  ONLY: ptable
  USE simulation_cell,                 ONLY: pbc
  USE string_utilities,                ONLY: str_search
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE topology_types,                  ONLY: topology_parameters_type

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "topology_util"


  PRIVATE
  PUBLIC :: find_boundary,&
            find_boundary2,&
            find_boundary3,&
            sort_bond,&
            sort_bend,&
            sort_dihe,&
            topology_connectivity_pack,&
            topology_coordinate_pack,&
            topology_generate_bend,&
            topology_generate_bond,&
            topology_generate_dihe,&
            topology_generate_molecule

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     topology_generate_molecule
!!
!!
!!   FUNCTION
!!     Use information from bond list to generate molecule. (ie clustering)
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE topology_generate_molecule(topology,globenv)


    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_molecule', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: first, iatom, iiatom, iii, &
                                                inum, istat, itype, iw, last, &
                                                natom, nbond, nnum, ntype
    INTEGER, POINTER                         :: checkme(:)
    LOGICAL                                  :: PM1, PM2, PM3, PM4

!# of molecular types
!# of molecules of each type
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

  IF(PM1)  WRITE(iw,*) "    Entering topology_generate_molecule"

  ntype = 1
  itype = 0
  nnum = 1
  inum = 0

  natom = topology%natoms
  nbond = SIZE(topology%bond_a)
 
  ALLOCATE(checkme(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_generate_molecule','checkme')
  ALLOCATE (topology%map_mol_typ(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_generate_molecule','map_mol_typ')
  ALLOCATE (topology%map_mol_num(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_generate_molecule','map_mol_num')

  !Zero the arrays
  checkme(:) = 0
  topology%map_mol_typ(:) = 0
  topology%map_mol_num(:) = 0

  !Parse the atoms list to see how many different molecule types there are
  ntype = 1
  topology%map_mol_typ(1)=1
  DO iatom=2,natom
    IF(topology%label_resname(iatom-1)==topology%label_resname(iatom)) THEN
      topology%map_mol_typ(iatom)=ntype
    ELSE
      ntype=ntype+1
      topology%map_mol_typ(iatom)=ntype
    END IF
  END DO

  DO itype=1,ntype
    !search for first atom in the type list and last to cut down time
    CALL find_boundary(topology%map_mol_typ,natom,first,last,itype)
    IF(PM2) WRITE(iw,*) "      itype ",itype,first,last

    inum = 0

    DO iatom=first,last
      checkme(:) = 0

      IF(topology%map_mol_num(iatom)==0) THEN
        inum=inum+1
        checkme(iatom)=1
        DO WHILE (SUM(checkme)>0)
          DO iiatom=first,last
            IF((checkme(iiatom)==1)) THEN
              checkme(iiatom)=0
              topology%map_mol_num(iiatom)=inum
              !-------------------------------------------------------------
              !-------------------------------------------------------------
              !Shortcut version comment out below is still buggy
              !CALL find_boundary(bond_a,nbond,a,b,iiatom)
              !WRITE(*,*) "Local find boundary",iiatom,a,b
              !DO iii=a,b
              !  checkme(bond_b(iii))=1
              !  WRITE(*,*) "iii = ",iii
              !END DO
              !-------------------------------------------------------------
              !-------------------------------------------------------------
              !Way slow version but can be improved with above
              DO iii=1,nbond
                IF((topology%bond_a(iii)==iiatom).AND.&
                  (topology%map_mol_num(topology%bond_b(iii))==0)) THEN 
                  !WRITE(*,*) "bond1 ",iiatom,topology%bond_b(iii)
                  checkme(topology%bond_b(iii))=1
                END IF
                IF((topology%bond_b(iii)==iiatom).AND.&
                  (topology%map_mol_num(topology%bond_a(iii))==0)) THEN
                  !WRITE(*,*) "bond2 ",iiatom,topology%bond_a(iii)
                  checkme(topology%bond_a(iii))=1
                END IF
              END DO
              !-------------------------------------------------------------
              !-------------------------------------------------------------
            END IF
          END DO
        END DO
      END IF

    END DO
  END DO

  DEALLOCATE(checkme,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_generate_molecule','checkme')

  IF(PM4) THEN
    DO iatom=1,natom
      WRITE(iw,*) "      iatom,map_mol_typ,map_mol_num ",iatom,&
                 topology%map_mol_typ(iatom),topology%map_mol_num(iatom)
    END DO
  END IF

  IF(PM1) WRITE(iw,*) "    Exiting  topology_generate_molecule"

  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("leaving "//routineN,globenv)
  END IF
END SUBROUTINE topology_generate_molecule

!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     find_boundary
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     find get the indices for the element that matches in the array
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE find_boundary(num_array,ntot,first,last,search)


    INTEGER, POINTER                         :: num_array(:)
    INTEGER                                  :: ntot, first, last, search

    INTEGER                                  :: i
    LOGICAL                                  :: found

!------------------------------------------------------------------------------

  found=.FALSE.
  first=0
  last=ntot

  DO i=1,ntot
    IF(num_array(i)==search) THEN
      IF(.NOT.found) THEN
        first = i
      END IF
      found=.TRUE.
    ELSE
      IF(found) THEN
        last = i-1
        EXIT                  !Exit this do loop???
      END IF
      found=.FALSE.
    END IF
  END DO

END SUBROUTINE find_boundary

!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     find_boundary2
!!
!!   FUNCTION
!!     find get the indices for the element that matches in the array
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE find_boundary2(num_array1,num_array2,ntot,first,last,search1,search2)


    INTEGER, POINTER                         :: num_array1(:), num_array2(:)
    INTEGER                                  :: ntot, first, last, search1, &
                                                search2

    INTEGER                                  :: i, tfirst, tlast
    LOGICAL                                  :: found

!------------------------------------------------------------------------------

  found=.FALSE.
  first=0
  last=ntot

  CALL find_boundary(num_array1,ntot,tfirst,tlast,search1)

  DO i=tfirst,tlast
    IF(num_array2(i)==search2) THEN
      IF(.NOT.found) THEN
        first = i
      END IF
      found=.TRUE.
    ELSE
      IF(found) THEN
        last = i-1
        EXIT                  !Exit this do loop???
      END IF
      found=.FALSE.
    END IF
  END DO

END SUBROUTINE find_boundary2

!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     find_boundary3
!!
!!   FUNCTION
!!     find get the indices for the element that matches in the array
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE find_boundary3(num_array1,num_array2,num_array3,ntot,first,last,&
                          search1,search2,search3)


    INTEGER, POINTER                         :: num_array1(:), num_array2(:), &
                                                num_array3(:)
    INTEGER                                  :: ntot, first, last, search1, &
                                                search2, search3

    INTEGER                                  :: i, tfirst, tlast
    LOGICAL                                  :: found

!------------------------------------------------------------------------------

  found=.FALSE.
  first=0
  last=ntot

  CALL find_boundary2(num_array1,num_array2,ntot,tfirst,tlast,search1,search2)

  DO i=tfirst,tlast
    IF(num_array3(i)==search3) THEN
      IF(.NOT.found) THEN
        first = i
      END IF
      found=.TRUE.
    ELSE
      IF(found) THEN
        last = i-1
        EXIT                  !Exit this do loop???
      END IF
      found=.FALSE.
    END IF
  END DO

END SUBROUTINE find_boundary3

!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     sort_bond
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     sort_the bonds lists
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE sort_bond(nbond,bond_a,bond_b)


    INTEGER                                  :: nbond
    INTEGER, POINTER                         :: bond_a(:), bond_b(:)

    INTEGER                                  :: first, i, j, k, last, search

!------------------------------------------------------------------------------

  DO i=1,nbond
    IF(bond_a(i)>bond_b(i)) THEN
      CALL swap(bond_a(i),bond_b(i))
    END IF
  END DO
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

END SUBROUTINE sort_bond

!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     sort_bend
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     sort_the bends lists
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE sort_bend(ntheta,theta_a,theta_b,theta_c)


    INTEGER                                  :: ntheta
    INTEGER, POINTER                         :: theta_a(:), theta_b(:), &
                                                theta_c(:)

    INTEGER                                  :: first, i, j, k, last, &
                                                search1, search2

!------------------------------------------------------------------------------

  DO i=1,ntheta
    IF(theta_a(i)>theta_c(i)) THEN
      CALL swap(theta_a(i),theta_c(i))
    END IF
  END DO
  DO i=1,ntheta
    DO j=1,ntheta
      IF(theta_a(i)<theta_a(j)) THEN
        CALL swap(theta_a(i),theta_a(j))
        CALL swap(theta_b(i),theta_b(j))
        CALL swap(theta_c(i),theta_c(j))
      END IF
    END DO
  END DO
  i=1
  DO WHILE (i<ntheta)
    search1 = theta_a(i)
    CALL find_boundary(theta_a,ntheta,first,last,search1)
    IF(.NOT.(first==last)) THEN
      DO j=first,last
        DO k=first,last
          IF(theta_b(j)<theta_b(k)) THEN
            CALL swap(theta_a(j),theta_a(k))
            CALL swap(theta_b(j),theta_b(k))
            CALL swap(theta_c(j),theta_c(k))
          END IF
        END DO
      END DO
      i=last+1
    ELSE
      i=i+1
    END IF
  END DO
  i=1
  DO WHILE (i<ntheta)
    search1 = theta_a(i)
    search2 = theta_b(i)
    CALL find_boundary2(theta_a,theta_b,ntheta,first,last,search1,search2)
    IF(.NOT.(first==last)) THEN
      DO j=first,last
        DO k=first,last
          IF(theta_c(j)<theta_c(k)) THEN
            CALL swap(theta_a(j),theta_a(k))
            CALL swap(theta_b(j),theta_b(k))
            CALL swap(theta_c(j),theta_c(k))
          END IF
        END DO
      END DO
      i=last+1
    ELSE
      i=i+1
    END IF
  END DO

END SUBROUTINE sort_bend

!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     sort_dihe
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     sort_the dihe lists
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE sort_dihe(nphi,phi_a,phi_b,phi_c,phi_d)


    INTEGER                                  :: nphi
    INTEGER, POINTER                         :: phi_a(:), phi_b(:), phi_c(:), &
                                                phi_d(:)

    INTEGER                                  :: first, i, j, k, last, &
                                                search1, search2, search3

!------------------------------------------------------------------------------
!Internal ordering in each dihedral

  DO i=1,nphi
    IF(phi_a(i)>phi_d(i)) THEN
      CALL swap(phi_a(i),phi_d(i))
      CALL swap(phi_b(i),phi_c(i))
    END IF
  END DO
  !Sort the first column
  DO i=1,nphi
    DO j=1,nphi
      IF(phi_a(i)<phi_a(j)) THEN
        CALL swap(phi_a(i),phi_a(j))
        CALL swap(phi_b(i),phi_b(j))
        CALL swap(phi_c(i),phi_c(j))
        CALL swap(phi_d(i),phi_d(j))
      END IF
    END DO
  END DO
  !Sort the second column
  i=1
  DO WHILE (i<nphi)
    search1 = phi_a(i)
    CALL find_boundary(phi_a,nphi,first,last,search1)
    IF(.NOT.(first==last)) THEN
      DO j=first,last
        DO k=first,last
          IF(phi_b(j)<phi_b(k)) THEN
            CALL swap(phi_a(j),phi_a(k))
            CALL swap(phi_b(j),phi_b(k))
            CALL swap(phi_c(j),phi_c(k))
            CALL swap(phi_d(j),phi_d(k))
          END IF
        END DO
      END DO
      i=last+1
    ELSE
      i=i+1
    END IF
  END DO
  !Sort the third column
  i=1
  DO WHILE (i<nphi)
    search1 = phi_a(i)
    search2 = phi_b(i)
    CALL find_boundary2(phi_a,phi_b,nphi,first,last,search1,search2)
    IF(.NOT.(first==last)) THEN
      DO j=first,last
        DO k=first,last
          IF(phi_c(j)<phi_c(k)) THEN
            CALL swap(phi_a(j),phi_a(k))
            CALL swap(phi_b(j),phi_b(k))
            CALL swap(phi_c(j),phi_c(k))
          END IF
        END DO
      END DO
      i=last+1
    ELSE
      i=i+1
    END IF
  END DO
  !Sort the fourth column
  i=1
  DO WHILE (i<nphi)
    search1 = phi_a(i)
    search2 = phi_b(i)
    search3 = phi_c(i)
    CALL find_boundary3(phi_a,phi_b,phi_c,nphi,first,last,&
                        search1,search2,search3)
    IF(.NOT.(first==last)) THEN
      DO j=first,last
        DO k=first,last
          IF(phi_d(j)<phi_d(k)) THEN
            CALL swap(phi_a(j),phi_a(k))
            CALL swap(phi_b(j),phi_b(k))
            CALL swap(phi_c(j),phi_c(k))
          END IF
        END DO
      END DO
      i=last+1
    ELSE
      i=i+1
    END IF
  END DO

END SUBROUTINE sort_dihe

!
!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     topology_generate_bond
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Use info from periodic table and assumptions to generate bonds
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE topology_generate_bond(topology,globenv)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_bond', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock = 1000

    CHARACTER(LEN=2), POINTER                :: atm_symbol(:)
    INTEGER                                  :: first, i, iatom, ibond, &
                                                ielem, istat, iw, j, jatom, &
                                                k, last, natom, nbond, nelem, &
                                                search
    INTEGER, POINTER                         :: tbond_a(:), tbond_b(:)
    LOGICAL                                  :: PM1, PM2, PM3, PM4
    REAL(dbl)                                :: dr(3), r2, rbond, rbond2, &
                                                rr(3)
    REAL(dbl), POINTER                       :: covalent_radius(:), &
                                                vdw_radius(:)

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

  IF(PM1) WRITE(iw,*) "    Entering topology_generate_bond"

  natom = topology%natoms

  !Get the temporary array that we need
  ALLOCATE(tbond_a(nblock),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','tbond_a',nblock)
  ALLOCATE(tbond_b(nblock),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','tbond_b',nblock)
  ALLOCATE(atm_symbol(natom),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','atm_symbol',natom)
  ALLOCATE(covalent_radius(natom),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','covalent_radius',natom)
  ALLOCATE(vdw_radius(natom),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','vdw_radius',natom)
!------------------------------------------------------------------------------
  DO iatom=1,natom
    atm_symbol(iatom) = topology%label_atmname(iatom)
  END DO
!------------------------------------------------------------------------------
  !Set up the arrays that we need for this routine
  tbond_a(:) = 0
  tbond_b(:) = 0

  DO iatom=1,natom
    nelem = UBOUND(ptable,1)
    DO ielem=1,nelem 
      IF(ptable(ielem)%symbol==atm_symbol(iatom)) THEN

        covalent_radius(iatom) = ptable(ielem)%covalent_radius
        CALL convert_to_cp2k_units("ANGSTROM",length=covalent_radius(iatom))

        vdw_radius(iatom) = ptable(ielem)%vdw_radius
        CALL convert_to_cp2k_units("ANGSTROM",length=vdw_radius(iatom))

        topology%atm_mass(iatom) = ptable(ielem)% amass
        !Need to figure out the correct charge!!!
        topology%atm_charge(iatom) = 0.0_dbl
        IF(PM2) WRITE(iw,'(7X,A,A5)') "iatom = ",ptable(ielem)%symbol
        !EXIT
      END IF
    END DO
  END DO
!------------------------------------------------------------------------------
  !Find bonds according to periodic table
  nbond = 0
  DO iatom=1,natom
    DO jatom=iatom+1,natom
      dr(1) = topology%r(1,iatom) - topology%r(1,jatom)
      dr(2) = topology%r(2,iatom) - topology%r(2,jatom)
      dr(3) = topology%r(3,iatom) - topology%r(3,jatom)

      rr = pbc ( dr, topology%cell )
      r2 = DOT_PRODUCT(rr,rr)

      IF(r2 <= 0.01_dbl) THEN
        CALL stop_program ("topology_generate_bond",&
                           "bond distance between atoms less then 0.1")
      END IF

      IF(topology%bondparm_type=="COVALENT") THEN
        rbond = covalent_radius(iatom)+covalent_radius(jatom)
      ELSE IF(topology%bondparm_type=="VDW") THEN
        rbond = MAX(vdw_radius(iatom),vdw_radius(jatom))
      ELSE
        CALL stop_program ("topology_generate_bond",&
                           "illegal bondparm_type")
      END IF
      rbond2 = rbond*rbond

      !Look for the right factor to use...
      rbond2 = rbond2*(topology%bondparm_factor)

      !Test the distance to the sum of the covalent radius
      IF(r2 < rbond2) THEN
        nbond=nbond+1
        IF(nbond > SIZE(tbond_a)) THEN
          CALL reallocate(tbond_a,1,nbond+nblock)
          CALL reallocate(tbond_b,1,nbond+nblock)
        END IF
        tbond_a(nbond)=iatom
        tbond_b(nbond)=jatom
        IF(PM2) WRITE(iw,'(7X,A,I5,I5,I5)') "bond ",nbond,iatom,jatom

        !Check to make sure no bonds are allowed between two different types
        !of molecules specified in the input file
        IF(topology%label_resname(iatom) /= topology%label_resname(jatom)) THEN
          CALL stop_program ("topology_generate_bond",&
                             "bonds between different molecule types??? ")
        END IF

      END IF
    END DO
  END DO
!------------------------------------------------------------------------------
  !Copy the stuff from temporary array into the more permanent one
  ALLOCATE(topology%bond_a(nbond),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','bond_a',nbond)
  ALLOCATE(topology%bond_b(nbond),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','bond_b',nbond)
  DO ibond=1,nbond
    topology%bond_a(ibond) = tbond_a(ibond)
    topology%bond_b(ibond) = tbond_b(ibond)
  END DO
!------------------------------------------------------------------------------
  !Sort out the bond list
  DO ibond=1,nbond
    IF(topology%bond_a(ibond)>topology%bond_b(ibond)) THEN 
      CALL swap(topology%bond_a(ibond),topology%bond_b(ibond))
    END IF
  END DO
  !Need to update this to quicksort for systems greater then 20000 atoms
  DO i=1,nbond
    DO j=1,nbond
      IF(topology%bond_a(i)<topology%bond_a(j)) THEN
        CALL swap(topology%bond_a(i),topology%bond_a(j))
        CALL swap(topology%bond_b(i),topology%bond_b(j))
      END IF
    END DO
  END DO
  i=1
  DO WHILE (i<nbond)
    search = topology%bond_a(i)
    CALL find_boundary(topology%bond_a,nbond,first,last,search)
    IF(.NOT.(first==last)) THEN
      DO j=first,last
        DO k=first,last
          IF(topology%bond_b(j)<topology%bond_b(k)) THEN
            CALL swap(topology%bond_a(j),topology%bond_a(k))
            CALL swap(topology%bond_b(j),topology%bond_b(k))
          END IF
        END DO
      END DO
      i=last+1
    ELSE
      i=i+1
    END IF
  END DO
!------------------------------------------------------------------------------
  !Sort out the duplicates... (Shouldn't have any but just to check)
  DO ibond=2,nbond
    IF((topology%bond_a(ibond-1)==topology%bond_a(ibond)).AND.&
       (topology%bond_b(ibond-1)==topology%bond_b(ibond))) THEN
      CALL stop_program ("topology_generate_bond",&
                         "duplicate bonds")
    END IF
  END DO

!------------------------------------------------------------------------------
  DEALLOCATE(tbond_a,STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','tbond_a')
  DEALLOCATE(tbond_b,STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','tbond_b')
  DEALLOCATE(atm_symbol,STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','atm_symbol')
  DEALLOCATE(covalent_radius,STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','covalent_radius')
  DEALLOCATE(vdw_radius,STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','vdw_radius')

  IF(PM1) WRITE(iw,*) "    Exiting  topology_generate_bond"

  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("leaving "//routineN,globenv)
  END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE topology_generate_bond


!
!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     topology_generate_bend
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Using a list of bonds, generate a list of bends
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE topology_generate_bend(topology,globenv)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_bend', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock = 1000

    INTEGER                                  :: countme, first, ibond, &
                                                iibond, istat, itheta, iw, &
                                                last, nbond, ntheta, search
    INTEGER, POINTER                         :: tmp_a(:), tmp_b(:), tmp_c(:)
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

  IF(PM1) WRITE(iw,*) "    Entering topology_generate_bend"

  nbond = SIZE(topology%bond_a)
  !WRITE(*,*) "nbond = ",nbond
  ntheta = 0
  !Stupid and inefficient but faster than alternative
  ALLOCATE(tmp_a(nblock),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bend','tmp_a',nblock)
  ALLOCATE(tmp_b(nblock),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bend','tmp_b',nblock)
  ALLOCATE(tmp_c(nblock),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bend','tmp_c',nblock)
!------------------------------------------------------------------------------
  ntheta=0
  DO ibond=1,nbond
    !Case A
    search = topology%bond_a(ibond)
    CALL find_boundary(topology%bond_a,nbond,first,last,search)
    IF (first /= 0) THEN
      DO iibond=first,last
        IF (topology%bond_b(ibond)/=topology%bond_b(iibond)) THEN
          ntheta = ntheta + 1
          IF(ntheta > SIZE(tmp_a)) THEN
            CALL reallocate(tmp_a,1,ntheta+nblock)
            CALL reallocate(tmp_b,1,ntheta+nblock)
            CALL reallocate(tmp_c,1,ntheta+nblock)
          END IF
          tmp_a(ntheta) = topology%bond_b(ibond)
          tmp_b(ntheta) = search
          tmp_c(ntheta) = topology%bond_b(iibond)
        END IF
      END DO
    END IF
    !Case B
    search = topology%bond_b(ibond)
    CALL find_boundary(topology%bond_a,nbond,first,last,search)
    IF (first /= 0) THEN
      DO iibond=first,last
        IF(topology%bond_a(ibond)/=topology%bond_b(iibond)) THEN
          ntheta = ntheta + 1
          IF(ntheta > SIZE(tmp_a)) THEN
            CALL reallocate(tmp_a,1,ntheta+nblock)
            CALL reallocate(tmp_b,1,ntheta+nblock)
            CALL reallocate(tmp_c,1,ntheta+nblock)
          END IF
          tmp_a(ntheta) = topology%bond_a(ibond)
          tmp_b(ntheta) = search
          tmp_c(ntheta) = topology%bond_b(iibond)
        END IF
      END DO
    END IF
    !Case C
    search = topology%bond_b(ibond)
    CALL find_boundary(topology%bond_b,nbond,first,last,search)
    IF (first /= 0) THEN
      DO iibond=first,last
        IF (topology%bond_a(ibond)/=topology%bond_a(iibond)) THEN
          ntheta = ntheta + 1
          IF(ntheta > SIZE(tmp_a)) THEN
            CALL reallocate(tmp_a,1,ntheta+nblock)
            CALL reallocate(tmp_b,1,ntheta+nblock)
            CALL reallocate(tmp_c,1,ntheta+nblock)
          END IF
          tmp_a(ntheta) = topology%bond_a(ibond)
          tmp_b(ntheta) = search
          tmp_c(ntheta) = topology%bond_a(iibond)
        END IF
      END DO
    END IF
  END DO
!------------------------------------------------------------------------------
  !Sort all the bends first
  CALL sort_bend(ntheta,tmp_a,tmp_b,tmp_c)

  ALLOCATE(topology%theta_a(ntheta),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bend','theta_a',ntheta)
  ALLOCATE(topology%theta_b(ntheta),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bend','theta_b',ntheta)
  ALLOCATE(topology%theta_c(ntheta),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bend','theta_c',ntheta)
  countme=1
  IF (ntheta > 0) THEN
    topology%theta_a(countme) = tmp_a(1)
    topology%theta_b(countme) = tmp_b(1)
    topology%theta_c(countme) = tmp_c(1)
    DO itheta=2,ntheta
      IF( .NOT.(tmp_a(itheta) == topology%theta_a(countme) .AND. &
          tmp_b(itheta) == topology%theta_b(countme) .AND. &
          tmp_c(itheta) == topology%theta_c(countme))) THEN
        countme = countme + 1
        topology%theta_a(countme) = tmp_a(itheta)
        topology%theta_b(countme) = tmp_b(itheta)
        topology%theta_c(countme) = tmp_c(itheta)
      END IF
    END DO
    ntheta = countme
    CALL REALLOCATE(topology%theta_a,1,ntheta)
    CALL REALLOCATE(topology%theta_b,1,ntheta)
    CALL REALLOCATE(topology%theta_c,1,ntheta)
  END IF
  

  !WRITE(*,*) "NTHETA = ",ntheta,SIZE(topology%theta_a)

!------------------------------------------------------------------------------
  DEALLOCATE(tmp_a)
  DEALLOCATE(tmp_b)
  DEALLOCATE(tmp_c)

  IF(PM1) WRITE(iw,*) "    Exiting  topology_generate_bend"

  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("leaving "//routineN,globenv)
  END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

END SUBROUTINE topology_generate_bend
 
 
 
!
!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     topology_generate_dihe
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Using a list of dihe, generate a list of dihe
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
!Use information from bond list to generate a list of dihedral angles
SUBROUTINE topology_generate_dihe(topology,globenv)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_generate_dihe', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iphi, iw, nphi, ntheta
    INTEGER, POINTER                         :: tmp_a(:), tmp_b(:), tmp_c(:), &
                                                tmp_d(:)
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

  IF(PM1) WRITE(iw,*) "    Entering topology_generate_dihe"

  ntheta = SIZE(topology%theta_a)
  nphi = 0
  !Stupid and inefficient but faster than alternative
  ALLOCATE(tmp_a(2*ntheta))
  ALLOCATE(tmp_b(2*ntheta))
  ALLOCATE(tmp_c(2*ntheta))
  ALLOCATE(tmp_d(2*ntheta))
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  ALLOCATE(topology%phi_a(nphi))
  ALLOCATE(topology%phi_b(nphi))
  ALLOCATE(topology%phi_c(nphi))
  ALLOCATE(topology%phi_d(nphi))
  DO iphi=1,nphi
    topology%phi_a(iphi) = tmp_a(iphi)
    topology%phi_b(iphi) = tmp_b(iphi)
    topology%phi_c(iphi) = tmp_c(iphi)
    topology%phi_d(iphi) = tmp_d(iphi)
  END DO
  nphi = SIZE(topology%phi_a)
!------------------------------------------------------------------------------
  DEALLOCATE(tmp_a)
  DEALLOCATE(tmp_b)
  DEALLOCATE(tmp_c)
  DEALLOCATE(tmp_d)


  IF(PM1) WRITE(iw,*) "    Exiting  topology_generate_dihe"

  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("leaving "//routineN,globenv)
  END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE topology_generate_dihe
 


!!*****
!******************************************************************************
!!****** topology/topology_connectivity_pack [1.0] *
!!
!!   NAME
!!     topology_connectivity_pack
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     topology connectivity pack
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE topology_connectivity_pack(molecule_kind_set,molecule_set,&
                                      topology,globenv)


    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_connectivity_pack', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: nhcopt
    INTEGER                                  :: atm_a, atm_b, atm_c, atm_d, &
                                                counter, first, i, istat, iw, &
                                                j, last, natom, nbond, nmol, &
                                                nphi, ntheta
    INTEGER, DIMENSION(:), POINTER           :: molecule_list
    INTEGER, POINTER                         :: checkme(:), first_list(:)
    LOGICAL                                  :: PM1, PM2, PM3, PM4
    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
    TYPE(bend_type), DIMENSION(:), POINTER   :: bend_list
    TYPE(bond_type), DIMENSION(:), POINTER   :: bond_list
    TYPE(local_molecule_type), POINTER       :: lmi
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(torsion_type), DIMENSION(:), &
      POINTER                                :: torsion_list

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

  IF(PM1) WRITE(iw,*) "  Entering topology_connectivity_pack"


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. Set the topology%[nmol_type,nmol]
  !-----------------------------------------------------------------------------
  natom  = topology%natoms
  topology%nmol_type = topology%map_mol_typ(natom)
  DO i=1,topology%nmol_type
    CALL find_boundary(topology%map_mol_typ,natom,first,last,i)
    topology%nmol=topology%nmol+topology%map_mol_num(last)
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. Allocate the molecule_kind_set
  !-----------------------------------------------------------------------------
  IF(topology%nmol_type <= 0) THEN
    CALL stop_program("topology_connectivity_pack","no molecule kind defined")
  ELSE
    NULLIFY(molecule_kind_set)
    CALL allocate_molecule_kind_set(molecule_kind_set,topology%nmol_type)
    IF(PM2) WRITE(iw,*) "    Allocated molecule_kind_set, Dimenstion of ",&
                  topology%nmol_type
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 3. Allocate the molecule_set
  !-----------------------------------------------------------------------------
  IF(topology%nmol <= 0) THEN
    CALL stop_program("topology_connectivity_pack","no molecule defined")
  ELSE
    NULLIFY(molecule_set)
    CALL allocate_molecule_set(molecule_set,topology%nmol)
    IF(PM2) WRITE(iw,*) "    Allocated molecule_set, dimenstion of ",&
                  topology%nmol
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 4. Set the molecule_kind%[kind_number,name,nhcopt]
  !-----------------------------------------------------------------------------
  counter=0
  nhcopt = topology % nhcopt
  DO i=1,topology%nmol_type
    natom  = topology%natoms
    molecule_kind => molecule_kind_set(i)
    CALL find_boundary(topology%map_mol_typ,natom,first,last,i)
    CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           kind_number=i,&
                           name=topology%label_resname(first),&
                           nhcopt=nhcopt)
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 5. Set the molecule_list for molecule_kind in molecule_kind_set
  !-----------------------------------------------------------------------------
  counter=0
  DO i=1,topology%nmol_type
    CALL find_boundary(topology%map_mol_typ,natom,first,last,i)
    nmol=topology%map_mol_num(last)
    ALLOCATE(molecule_list(nmol),STAT=istat)
    IF(istat/=0) CALL stop_memory('coordinate_pack','molecule_list',nmol)
    IF(PM2) WRITE(iw,'(5X,A,I5,I5)') "ikind,nmol",i,nmol
    DO j=1,nmol
      counter=counter+1
      molecule_list(j)=counter
    END DO
    molecule_kind => molecule_kind_set(i)
    CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           molecule_list=molecule_list)
    IF(PM2) WRITE(iw,*) "      molecule_list",molecule_list(:)
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 6. Set the molecule_set(imol)%first_atom via set_molecule_set
  !-----------------------------------------------------------------------------
  ALLOCATE(first_list(topology%nmol),STAT=istat)
  IF(istat/=0) CALL stop_memory('coordinate_pack','frist_list',topology%nmol)
  counter=0
  DO i=1,topology%nmol_type
    CALL find_boundary(topology%map_mol_typ,natom,first,last,i)
    nmol=topology%map_mol_num(last)
    DO j=1,nmol
      counter=counter+1
      CALL find_boundary2(topology%map_mol_typ,topology%map_mol_num,&
                          natom,first,last,i,j)
      first_list(counter)=first
      IF(PM2) WRITE(iw,'(5X,A,I5,I5,I5)') "first_list=ikind,nmol,first_atom",&
                    i,j,first
    END DO
  END DO
  CALL set_molecule_set(molecule_set,first_list)
  DEALLOCATE(first_list,STAT=istat)
  IF(istat/=0) CALL stop_memory('coordinate_pack','frist_list')


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 7. Set and NULLIFY the molecule_set(imol)%lmi via set_molecule_set
  !-----------------------------------------------------------------------------
  DO i=1,topology%nmol
    molecule => molecule_set(i)
    NULLIFY(lmi)
    ALLOCATE(lmi,STAT=istat)
    IF(istat/=0) CALL stop_memory('connectivity_pack','lmi',1)
    lmi%nstates=0
    NULLIFY(lmi%states)
    CALL set_molecule(molecule,lmi=lmi)
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 8. Set the atom_list for molecule_kind in molecule_kind_set (PART 1)
  !-----------------------------------------------------------------------------
  counter=0
  DO i=1,topology%nmol_type
    molecule_kind => molecule_kind_set(i)
    CALL find_boundary(topology%map_mol_typ,topology%natoms,first,last,i)
    natom = (last-first+1)/ topology%map_mol_num(last)
    ALLOCATE(atom_list(natom),STAT=istat)
    IF(istat/=0) CALL stop_memory('connectivity_pack','atom_list',natom)
    DO j=1,natom
      !Atomic kind information will be filled in (PART 2)
      NULLIFY(atom_list(j)%atomic_kind)
      atom_list(j)%name=topology%label_atmname(j+first-1)
      IF(PM2) WRITE(iw,'(5X,A,I5,I5,1X,A5)') "atom_list ",i,j,atom_list(j)%name
    END DO
    CALL set_molecule_kind(molecule_kind=molecule_kind,atom_list=atom_list)
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 9. Set the molecule_set(imol)%molecule_kind via set_molecule
  !-----------------------------------------------------------------------------
  counter=0
  DO i=1,topology%nmol_type
    molecule_kind => molecule_kind_set(i)
    CALL find_boundary(topology%map_mol_typ,topology%natoms,first,last,i)
    nmol=topology%map_mol_num(last)
    DO j=1,nmol
      counter=counter+1
      molecule => molecule_set(counter)
      CALL set_molecule(molecule,molecule_kind=molecule_kind)
    END DO
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 10. Set the molecule_kind%[nbond,bond_list] via set_molecule_kind
  !-----------------------------------------------------------------------------
  DO i=1,topology%nmol_type
    molecule_kind => molecule_kind_set(i)
    nbond  = SIZE(topology%bond_a)
    CALL find_boundary(topology%map_mol_typ,topology%natoms,first,last,i)
    ALLOCATE(checkme(nbond),STAT=istat)
    IF(istat/=0) CALL stop_memory('connectivity_pack','checkme',nbond)
    checkme(:)=0
    DO j=1,nbond
      atm_a = topology%bond_a(j)
      atm_b = topology%bond_b(j)
      IF( (topology%map_mol_typ(atm_a)==i).AND.&
          (topology%map_mol_num(atm_a)==1).AND.&
          (topology%map_mol_typ(atm_b)==i).AND.&
          (topology%map_mol_num(atm_b)==1) ) THEN
        checkme(j) = 1
      END IF
    END DO
    IF(PM2) WRITE(iw,*) "    Total number bonds for molecule type ",SUM(checkme)
    ALLOCATE(bond_list(SUM(checkme)),STAT=istat)
    IF(istat/=0) CALL stop_memory('connectivity_pack','bond_list',SUM(checkme))
    counter=0
    DO j=1,SUM(checkme)
      IF (checkme(j)==1) THEN
        counter=counter+1
        bond_list(counter)%a=topology%bond_a(j)-first+1
        bond_list(counter)%b=topology%bond_b(j)-first+1
        !point this to the right bond_kind_type if using force field
        NULLIFY(bond_list(j)%bond_kind)
        IF(PM2) THEN
          WRITE(iw,'(7X,A,I3,1X,A,I5,I5,1X,A,I5,I5)') "molecule_kind",i,"bond",&
                      topology%bond_a(j),&
                      topology%bond_b(j),&
                      "offset number at",&
                      topology%bond_a(j)-first+1,&
                      topology%bond_b(j)-first+1
        END IF
      END IF
    END DO
    CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           nbond=SUM(checkme),bond_list=bond_list)
    DEALLOCATE(checkme,STAT=istat)
    IF (istat /= 0) CALL stop_memory ('connectivity_pack','checkme')
  END DO
  counter=0
  DO i=1,topology%nmol_type
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           nmolecule=j,nbond=nbond)
    counter=counter+nbond*j
  END DO
  IF(counter/=SIZE(topology%bond_a)) THEN
    IF(PM2) WRITE(iw,*) "    checksum vs nbond ",counter,SiZE(topology%bond_a)
    CALL stop_program ( "connectivity_pack", "nbond counting error")
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 11. Set the molecule_kind%[nbend,bend_list] via set_molecule_kind
  !-----------------------------------------------------------------------------
  DO i=1,topology%nmol_type
    molecule_kind => molecule_kind_set(i)
    ntheta  = SIZE(topology%theta_a)
    CALL find_boundary(topology%map_mol_typ,topology%natoms,first,last,i)
    ALLOCATE(checkme(ntheta),STAT=istat)
    IF(istat/=0) CALL stop_memory('connectivity_pack','checkme',ntheta)
    checkme(:)=0
    DO j=1,ntheta
      atm_a = topology%theta_a(j)
      atm_b = topology%theta_b(j)
      atm_c = topology%theta_c(j)
      IF( (topology%map_mol_typ(atm_a)==i).AND.&
          (topology%map_mol_num(atm_a)==1).AND.&
          (topology%map_mol_typ(atm_b)==i).AND.&
          (topology%map_mol_num(atm_b)==1).AND.&
          (topology%map_mol_typ(atm_c)==i).AND.&
          (topology%map_mol_num(atm_c)==1) ) THEN
        checkme(j) = 1
      END IF
    END DO
    IF(PM2) WRITE(iw,*) "    Total number bends for molecule type ",SUM(checkme)
    ALLOCATE(bend_list(SUM(checkme)),STAT=istat)
    IF(istat/=0) CALL stop_memory('connectivity_pack','bend_list',SUM(checkme))
    counter=0
    DO j=1,SUM(checkme)
      IF (checkme(j)==1) THEN
        counter=counter+1
        bend_list(counter)%a=topology%theta_a(j)-first+1
        bend_list(counter)%b=topology%theta_b(j)-first+1
        bend_list(counter)%c=topology%theta_c(j)-first+1
        !point this to the right bond_kind_type if using force field
        NULLIFY(bend_list(j)%bend_kind)
        IF(PM2) THEN
          WRITE(iw,'(7X,A,I3,1X,A,I5,I5,I5,1X,A,I5,I5,I5)') &
                      "molecule_kind",i,"bend",&
                      topology%theta_a(j),&
                      topology%theta_b(j),&
                      topology%theta_c(j),&
                      "offset number at",&
                      topology%theta_a(j)-first+1,&
                      topology%theta_b(j)-first+1,&
                      topology%theta_c(j)-first+1
        END IF
      END IF
    END DO
    CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           nbend=SUM(checkme),bend_list=bend_list)
    DEALLOCATE(checkme,STAT=istat)
    IF (istat /= 0) CALL stop_memory ('connectivity_pack','checkme')
  END DO
  counter=0
  DO i=1,topology%nmol_type
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           nmolecule=j,nbend=ntheta)
    counter=counter+ntheta*j
  END DO
  IF(counter/=SIZE(topology%theta_a)) THEN
    IF(PM2) WRITE(iw,*) "    checksum vs nbend ",counter,SiZE(topology%theta_a)
    CALL stop_program ( "connectivity_pack", "nbend counting error")
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 12. Set the molecule_kind%[ntorsion,torsion_list] via set_molecule_kind
  !-----------------------------------------------------------------------------
  DO i=1,topology%nmol_type
    molecule_kind => molecule_kind_set(i)
    nphi  = SIZE(topology%phi_a)
    CALL find_boundary(topology%map_mol_typ,topology%natoms,first,last,i)
    ALLOCATE(checkme(nphi),STAT=istat)
    IF(istat/=0) CALL stop_memory('connectivity_pack','checkme',nphi)
    checkme(:)=0
    DO j=1,nphi
      atm_a = topology%phi_a(j)
      atm_b = topology%phi_b(j)
      atm_c = topology%phi_c(j)
      atm_d = topology%phi_d(j)
      IF( (topology%map_mol_typ(atm_a)==i).AND.&
          (topology%map_mol_num(atm_a)==1).AND.&
          (topology%map_mol_typ(atm_b)==i).AND.&
          (topology%map_mol_num(atm_b)==1).AND.&
          (topology%map_mol_typ(atm_c)==i).AND.&
          (topology%map_mol_num(atm_c)==1).AND.&
          (topology%map_mol_typ(atm_d)==i).AND.&
          (topology%map_mol_num(atm_d)==1) ) THEN
        checkme(j) = 1
      END IF
    END DO
    IF(PM2) WRITE(iw,*) "    Total torsions for molecule type ",SUM(checkme)
    ALLOCATE(torsion_list(SUM(checkme)),STAT=istat)
    IF(istat/=0) &
      CALL stop_memory('connectivity_pack','torsion_list',SUM(checkme))
    counter=0
    DO j=1,SUM(checkme)
      IF (checkme(j)==1) THEN
        counter=counter+1
        torsion_list(counter)%a=topology%phi_a(j)-first+1
        torsion_list(counter)%b=topology%phi_b(j)-first+1
        torsion_list(counter)%c=topology%phi_c(j)-first+1
        torsion_list(counter)%d=topology%phi_d(j)-first+1
        !point this to the right bond_kind_type if using force field
        NULLIFY(torsion_list(j)%torsion_kind)
        IF(PM2) THEN
          WRITE(iw,'(7X,A,I3,1X,A,I5,I5,I5,I5,1X,A,I5,I5,I5,I5)') &
                      "molecule_kind",i,"torsion",&
                      topology%phi_a(j),&
                      topology%phi_b(j),&
                      topology%phi_c(j),&
                      topology%phi_d(j),&
                      "offset number at",&
                      topology%phi_a(j)-first+1,&
                      topology%phi_b(j)-first+1,&
                      topology%phi_c(j)-first+1,&
                      topology%phi_d(j)-first+1
        END IF
      END IF
    END DO
    CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           ntorsion=SUM(checkme),torsion_list=torsion_list)
    DEALLOCATE(checkme,STAT=istat)
    IF (istat /= 0) CALL stop_memory ('connectivity_pack','checkme')
  END DO
  counter=0
  DO i=1,topology%nmol_type
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           nmolecule=j,ntorsion=nphi)
    counter=counter+nphi*j
  END DO
  IF(counter/=SIZE(topology%phi_a)) THEN
    IF(PM2) WRITE(iw,*) "    checksum vs nphi ",counter,SiZE(topology%phi_a)
    CALL stop_program ( "connectivity_pack", "nphi counting error")
  END IF


    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Stuff the ndc info
    !molecule_set(i)%molpar%ndcon = 0
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Stuff the n3x3c info
    !molecule_set(i)%molpar%n3x3con = 0
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------


  IF(PM1) WRITE(iw,*) "  Exiting  topology_connectivity_pack"

  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("leaving "//routineN,globenv)
  END IF

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
END SUBROUTINE topology_connectivity_pack




!******************************************************************************
!!****** topology_util/topology_coordinate_pack [1.0] *
!!
!!   NAME
!!     topology_coordinate_pack
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Take info readin from different file format and stuff it into
!!     compatible data structure in cp2k
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*** *************************************************************************

SUBROUTINE topology_coordinate_pack (particle_set,atomic_kind_set,&
                              molecule_kind_set,molecule_set,&
                              topology,globenv)


    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_coordinate_pack', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      POINTER                                :: work(:)
    INTEGER                                  :: counter, first, i, ikind, &
                                                istat, iw, j, k, last, natom
    INTEGER, DIMENSION(:), POINTER           :: atomlist, kind_of, &
                                                natom_of_kind
    LOGICAL                                  :: PM1, PM2, PM3, PM4
    REAL(dbl), DIMENSION(:), POINTER         :: charge, mass
    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(elp_potential_type), POINTER        :: elp_potential
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

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

  IF(PM1) WRITE(iw,*) "  Entering topology_coordinate_pack"


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. Determine topology%[natom_type,atom_names] and save mass(natom_type) 
  !-----------------------------------------------------------------------------
  counter = 0
  ALLOCATE(work(topology%natoms),mass(topology%natoms),STAT=istat)
  IF (istat/=0) CALL stop_memory ('coordinate_pack','work,mass',topology%natoms)
  work(:)=""
  DO i = 1, topology%nmol_type
    CALL find_boundary(topology%map_mol_typ,topology%natoms,first,last,i)
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           natom=natom,atom_list=atom_list)
    DO j = 1, natom
      IF(str_search(work,counter,atom_list(j)%name)==0)THEN
         counter = counter + 1
         work(counter) = atom_list(j)%name
         mass(counter) = topology%atm_mass(first+j-1)
        IF(PM2) WRITE(iw,'(7X,A,1X,A5,F10.5)') &
                "NEW ATOMIC KIND",work(counter),mass(counter)
      END IF
    END DO
  END DO
  topology%natom_type = counter
  ALLOCATE(topology%atom_names(topology%natom_type),STAT=istat)
  IF(istat/=0) CALL stop_memory ('coordinate_pack','atom_names')
  topology%atom_names(1:counter)=work(1:counter)
  DEALLOCATE(work,STAT=istat)
  IF (istat/=0) CALL stop_memory ('coordinate_pack','work')
  CALL reallocate(mass,1,counter)
  IF(PM2) &
    WRITE(iw,'(5X,A,I3)') "Total Number of Atomic Kinds = ",topology%natom_type
  

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. Allocate the data structure for the atomic kind information 
  !-----------------------------------------------------------------------------
  NULLIFY(atomic_kind_set)
  CALL allocate_atomic_kind_set(atomic_kind_set,topology%natom_type)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 3.  Allocate the data structure for the atomic information
  !-----------------------------------------------------------------------------
  NULLIFY(particle_set)
  CALL allocate_particle_set(particle_set,topology%natoms)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 4. Set the atomic_kind_set(ikind)%[name,kind_number,mass]
  !-----------------------------------------------------------------------------
  DO i=1,topology%natom_type
    atomic_kind => atomic_kind_set(i)
    CALL set_atomic_kind(atomic_kind=atomic_kind,kind_number=i,&
                         name=topology%atom_names(i),mass=mass(i)) 
  END DO
  DEALLOCATE(mass,STAT=istat)
  IF (istat/=0) CALL stop_memory ('coordinate_pack','mass')


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 5. Determine number of atom of each kind (ie natom_of_kind and kind_of)
  !-----------------------------------------------------------------------------
  ALLOCATE(kind_of(topology%natoms),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','kind_of')
  ALLOCATE(natom_of_kind(topology%natom_type),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','natom_of_kind')
  kind_of(:) = 0
  natom_of_kind(:) = 0
  DO i=1,topology%natom_type
    DO j=1,topology%natoms
      IF(topology%atom_names(i)==topology%label_atmname(j)) THEN
        natom_of_kind(i)=natom_of_kind(i)+1
        IF(kind_of(j)==0) kind_of(j)=i
      END IF
    END DO
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 6. Set the atom_kind_set(ikind)%[natom,atom_list]
  !-----------------------------------------------------------------------------
  DO i=1,topology%natom_type
    atomic_kind => atomic_kind_set(i)
    NULLIFY (atomlist)
    ALLOCATE(atomlist(natom_of_kind(i)),STAT=istat)
    IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','atomlist')
    counter=0
    DO j=1,topology%natoms
      IF(kind_of(j)==i) THEN
        counter=counter+1
        atomlist(counter)=j
      END IF
    END DO
    IF(PM2) &  
      WRITE(iw,*) "      Atomic kind",i,"contain particle",atomlist(:)
    CALL set_atomic_kind(atomic_kind=atomic_kind,&
                         natom=natom_of_kind(i),atom_list=atomlist)
    DEALLOCATE(atomlist,STAT=istat)
    IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','atomlist')
  END DO
  DEALLOCATE(natom_of_kind,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','natom_of_kind')


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 7. Fill in some more infor for particle_set
  !-----------------------------------------------------------------------------
  DO i=1,topology%natoms
    ikind = kind_of(i)
    particle_set(i)%atomic_kind => atomic_kind_set(ikind)
    particle_set(i)%iatom = i
    particle_set(i)%r(:) = topology%r(:,i)
    NULLIFY(particle_set(i)%lpi)
  END DO
  DEALLOCATE(kind_of,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','kind_of')


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 8. Set the atom_list for molecule_kind in molecule_kind_set (PART 2)
  !-----------------------------------------------------------------------------
    DO i=1,topology%nmol_type
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           natom=natom,atom_list=atom_list)
    DO j=1,natom
      DO k=1,topology%natom_type
        IF(atom_list(j)%name == topology%atom_names(k)) THEN
          atom_list(j)%atomic_kind => atomic_kind_set(k)
        END IF
      END DO
    END DO
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 9. ????????????????????????????????????????????
  !-----------------------------------------------------------------------------
  ALLOCATE(charge(topology%natom_type),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','charge')
  charge(:)=0.0_dbl
  IF (globenv%program_name == "FIST") THEN
    CALL  init_atomic_kind_set(atomic_kind_set,globenv)
    DO ikind=1,topology%natom_type
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                        elp_potential=elp_potential) 
      CALL set_potential(potential=elp_potential,&
                        qeff=charge(ikind))
    END DO
  END IF
  DEALLOCATE (charge,STAT=istat)
  IF (istat /= 0) THEN
    CALL stop_memory('topology_coordinate_pack',&
              'release work storage')
  END IF

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  IF(PM1) WRITE(iw,*) "  Exiting  topology_coordinate_pack"
  CALL write_checkpoint_information("leaving "//routineN,globenv)
  !CALL stop_program("topology_coordinate_pack","checkpoint end")

END SUBROUTINE topology_coordinate_pack


END MODULE topology_util
