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
!!     Control for reading in different topologies
!!
!!   AUTHOR
!!     IKUO 08.01.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE topology_util

  USE atomic_kind_types,               ONLY: allocate_atomic_kind_set,&
                                             atomic_kind_type,&
                                             get_atomic_kind,&
                                             init_atomic_kind_set,&
                                             set_atomic_kind,&
                                             write_atomic_kind
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE external_potential_types,        ONLY: epc_potential_type,&
                                             set_potential
  USE global_types,                    ONLY: global_environment_type
  USE input_types,                     ONLY: setup_parameters_type
  USE kinds,                           ONLY: dbl
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_bcast
  USE molecule_types_new,              ONLY: allocate_molecule_set,&
                                             get_molecule,&
                                             local_molecule_type,&
                                             molecule_type,&
                                             set_molecule,&
                                             set_molecule_set
  USE molecule_kind_types,             ONLY: allocate_molecule_kind_set,&
                                             atom_type,&
                                             bond_type,&
                                             bend_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type,&
                                             read_molecule_kind_set,&
                                             set_molecule_kind,&
                                             torsion_type,&
                                             write_molecule_kind_set
  USE nrutil,                          ONLY: swap
  USE particle_types,                  ONLY: allocate_particle_set,&
                                             particle_type,&
                                             write_particle_coordinates
  USE periodic_table,                  ONLY: ptable
  USE simulation_cell,                 ONLY: cell_type,&
                                             write_cell,&
                                             pbc
  USE string_utilities,                ONLY: integer_to_string,&
                                             make_tuple,&
                                             str_search,&
                                             string_to_integer,&
                                             uppercase,&
                                             xstring
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE topology_types,                  ONLY: topology_parameters_type
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
!!   SYNOPSIS
!!     Subroutine topology_generate_molecule(natom, atm_mol_name, nbond,&
!!         bond_a, bond_b, map_mol_typ, map_mol_num)
!!       Implicit None
!!       Integer:: natom, nbond
!!       Character(Len=20), Pointer:: atm_mol_name(:)
!!       Integer, Pointer:: bond_a(:)
!!       Integer, Pointer:: bond_b(:)
!!       Integer, Pointer:: map_mol_typ(:)
!!       Integer, Pointer:: map_mol_num(:)
!!     End Subroutine topology_generate_molecule
!!
!!   FUNCTION
!!     Use information from bond list to generate molecule. (ie clustering)
!!
!!   AUTHOR
!!     IKUO 08.01.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE topology_generate_molecule(topology)

  IMPLICIT NONE

! Arguments
  TYPE ( topology_parameters_type ), INTENT (INOUT)   :: topology

! Locals
  INTEGER :: natom,nbond
  INTEGER :: iatom,iiatom,ibond
  INTEGER :: first,last
  INTEGER :: itype,ntype       !# of molecular types
  INTEGER :: inum,nnum         !# of molecules of each type
  INTEGER, POINTER :: checkme(:)
  INTEGER :: a,b,iii
  INTEGER :: istat
  CHARACTER (LEN=20) :: str1,str2

!------------------------------------------------------------------------------
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
                IF (topology%bond_a(iii)==iiatom) THEN
                  checkme(topology%bond_b(iii))=1
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

END SUBROUTINE topology_generate_molecule

!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     find_boundary
!!
!!   SYNOPSIS
!!     Subroutine find_boundary(num_array, ntot, first, last, search)
!!       Implicit None
!!       Integer, Pointer:: num_array(:)
!!       Integer:: ntot, first, last, search
!!     End Subroutine find_boundary
!!
!!   FUNCTION
!!     find get the indices for the element that matches in the array
!!
!!   AUTHOR
!!     IKUO 08.01.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE find_boundary(num_array,ntot,first,last,search)

  IMPLICIT NONE

! Arguments
  INTEGER, POINTER :: num_array(:)
  INTEGER :: ntot,first,last,search

! Locals
  LOGICAL :: found
  INTEGER :: i
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
!!     IKUO 08.17.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE find_boundary2(num_array1,num_array2,ntot,first,last,search1,search2)

  IMPLICIT NONE

! Arguments
  INTEGER, POINTER :: num_array1(:),num_array2(:)
  INTEGER :: ntot,first,last,search1,search2
  INTEGER :: tfirst,tlast

! Locals
  LOGICAL :: found
  INTEGER :: i
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
!!     IKUO 08.17.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE find_boundary3(num_array1,num_array2,num_array3,ntot,first,last,&
                          search1,search2,search3)

  IMPLICIT NONE

! Arguments
  INTEGER, POINTER :: num_array1(:),num_array2(:),num_array3(:)
  INTEGER :: ntot,first,last,search1,search2,search3
  INTEGER :: tfirst,tlast

! Locals
  LOGICAL :: found
  INTEGER :: i
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
!!     Subroutine sort_bond(nbond, bond_a, bond_b)
!!       Implicit None
!!       Integer:: nbond
!!       Integer, Pointer:: bond_a(:), bond_b(:)
!!     End Subroutine sort_bond
!!
!!   FUNCTION
!!     sort_the bonds lists
!!
!!   AUTHOR
!!     IKUO 08.17.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE sort_bond(nbond,bond_a,bond_b)

  IMPLICIT NONE

! Arguments
  INTEGER :: nbond
  INTEGER, POINTER :: bond_a(:),bond_b(:)

! Locals
  INTEGER :: i,j,k
  INTEGER :: first,last,search

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
!!     Subroutine sort_bend(ntheta, theta_a, theta_b, theta_c)
!!       Implicit None
!!       Integer:: ntheta
!!       Integer, Pointer:: theta_a(:), theta_b(:), theta_c(:)
!!     End Subroutine sort_bend
!!
!!   FUNCTION
!!     sort_the bends lists
!!
!!   AUTHOR
!!     IKUO 08.17.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE sort_bend(ntheta,theta_a,theta_b,theta_c)

  IMPLICIT NONE

! Arguments
  INTEGER :: ntheta
  INTEGER, POINTER :: theta_a(:),theta_b(:),theta_c(:)

! Locals
  INTEGER :: i,j,k
  INTEGER :: first,last,search1,search2

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
!!     Subroutine sort_dihe(nphi, phi_a, phi_b, phi_c, phi_d)
!!       Implicit None
!!       Integer:: nphi
!!       Integer, Pointer:: phi_a(:), phi_b(:), phi_c(:), phi_d(:)
!!     End Subroutine sort_dihe
!!
!!   FUNCTION
!!     sort_the dihe lists
!!
!!   AUTHOR
!!     IKUO 08.17.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE sort_dihe(nphi,phi_a,phi_b,phi_c,phi_d)

  IMPLICIT NONE

! Arguments
  INTEGER :: nphi
  INTEGER, POINTER :: phi_a(:),phi_b(:),phi_c(:),phi_d(:)

! Locals
  INTEGER :: i,j,k
  INTEGER :: first,last,search1,search2,search3

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
!!     Subroutine topology_generate_bond(cell, r, label_atmname, bond_a,&
!!         bond_b, atm_mass, atm_charge)
!!       Implicit None
!!       Type(cell_type), Pointer:: cell
!!       Real(Kind=dbl), Dimension(:,:), Pointer:: r
!!       Character(Len=20), Pointer:: label_atmname(:)
!!       Integer, Pointer:: bond_a(:), bond_b(:)
!!       Real(Kind=dbl), Pointer:: atm_charge(:)
!!       Real(Kind=dbl), Pointer:: atm_mass(:)
!!     End Subroutine topology_generate_bond
!!
!!   FUNCTION
!!     Use info from periodic table and assumptions to generate bonds
!!
!!   AUTHOR
!!     IKUO 08.01.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
SUBROUTINE topology_generate_bond(topology)
  IMPLICIT NONE
!Arguments
  TYPE ( topology_parameters_type ), INTENT (INOUT)   :: topology

!Locals
  INTEGER :: i,j,k
  INTEGER :: natom
  INTEGER :: nbond
  INTEGER :: iatom,jatom,ibond
  INTEGER :: first,last,search
  INTEGER :: ielem,nelem
  INTEGER, POINTER :: tbond_a(:),tbond_b(:)
  REAL(dbl) :: dr(3),rr(3),r2,rbond,rbond2
  REAL(dbl), POINTER :: pos_x(:),pos_y(:),pos_z(:)
  REAL(dbl), POINTER :: covalent_radius(:)
  CHARACTER (LEN = 2),POINTER :: atm_symbol(:)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  natom = topology%natoms

  WRITE(*,*) "natom = ",natom

  !Get the temporary array that we need
  ALLOCATE(tbond_a(2*natom))
  ALLOCATE(tbond_b(2*natom))
  ALLOCATE(pos_x(natom))
  ALLOCATE(pos_y(natom))
  ALLOCATE(pos_z(natom))
  ALLOCATE(atm_symbol(natom))
  ALLOCATE(covalent_radius(natom))
  ALLOCATE(topology%atm_mass(natom))
  ALLOCATE(topology%atm_charge(natom))
!------------------------------------------------------------------------------
  DO iatom=1,natom
    atm_symbol(iatom) = topology%label_atmname(iatom)
  END DO
!------------------------------------------------------------------------------
  !Set up the arrays that we need for this routine
  tbond_a(:) = 0
  tbond_b(:) = 0

  DO iatom=1,natom
    pos_x(iatom) = topology%r(1,iatom)
    pos_y(iatom) = topology%r(2,iatom)
    pos_z(iatom) = topology%r(3,iatom)
    
    nelem = SIZE(ptable)
    DO ielem=1,nelem 
      IF(ptable(ielem)%symbol==atm_symbol(iatom)) THEN
        covalent_radius(iatom) = ptable(ielem)%covalent_radius
        CALL convert_to_cp2k_units("ANGSTROM",length=covalent_radius(iatom))
        topology%atm_mass(iatom) = ptable(ielem)% amass
        !Need to figure out the correct charge!!!
        topology%atm_charge(iatom) = 0.0
        EXIT
      END IF
    END DO
  END DO
!------------------------------------------------------------------------------
  !Find bonds according to periodic table
  nbond = 0
  DO iatom=1,natom
    DO jatom=iatom+1,natom
      dr(1) = pos_x(iatom) - pos_x(jatom)
      dr(2) = pos_y(iatom) - pos_y(jatom)
      dr(3) = pos_z(iatom) - pos_z(jatom)
      rr = pbc ( dr, topology%cell )
      r2 = DOT_PRODUCT(rr,rr)
      rbond = covalent_radius(iatom)+covalent_radius(jatom)
      rbond2 = rbond*rbond
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !Look for the right factor to use...
      rbond2 = rbond2*1.1
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !Test the distance to the sum of the covalent radius
      IF(r2 < rbond2) THEN
        nbond=nbond+1
        tbond_a(nbond)=iatom
        tbond_b(nbond)=jatom
        !WRITE(*,*) "bond ",nbond,iatom,jatom
      END IF
    END DO
  END DO
!------------------------------------------------------------------------------
  !Copy the stuff from temporary array into the more permanent one
  ALLOCATE(topology%bond_a(nbond))
  ALLOCATE(topology%bond_b(nbond))
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
  DEALLOCATE(tbond_a)
  DEALLOCATE(tbond_b)
  DEALLOCATE(pos_x)
  DEALLOCATE(pos_y)
  DEALLOCATE(pos_z)
  DEALLOCATE(atm_symbol)
  DEALLOCATE(covalent_radius)
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
!!     Subroutine topology_generate_bend(bond_a, bond_b, theta_a, theta_b,&
!!         theta_c)
!!       Implicit None
!!       Integer, Pointer:: bond_a(:), bond_b(:)
!!       Integer, Pointer:: theta_a(:), theta_b(:), theta_c(:)
!!     End Subroutine topology_generate_bend
!!
!!   FUNCTION
!!     Using a list of bonds, generate a list of bends
!!
!!   AUTHOR
!!     IKUO 08.01.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE topology_generate_bend(topology)
  IMPLICIT NONE

!Arguments
  TYPE ( topology_parameters_type ), INTENT (INOUT)   :: topology

!Locals
  INTEGER :: ntheta
  INTEGER :: first,last,search,countme
  INTEGER :: nbond
  INTEGER :: ibond,iibond,itheta
  INTEGER, POINTER :: tmp_a(:),tmp_b(:),tmp_c(:)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  nbond = SIZE(topology%bond_a)
  !WRITE(*,*) "nbond = ",nbond
  ntheta = 0
  !Stupid and inefficient but faster than alternative
  ALLOCATE(tmp_a(2*nbond))
  ALLOCATE(tmp_b(2*nbond))
  ALLOCATE(tmp_c(2*nbond))
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

  ALLOCATE(topology%theta_a(ntheta))
  ALLOCATE(topology%theta_b(ntheta))
  ALLOCATE(topology%theta_c(ntheta))
  countme=1
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

  !WRITE(*,*) "NTHETA = ",ntheta

!------------------------------------------------------------------------------
  DEALLOCATE(tmp_a)
  DEALLOCATE(tmp_b)
  DEALLOCATE(tmp_c)
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
!!     Subroutine topology_generate_dihe(theta_a, theta_b, theta_c, phi_a,&
!!         phi_b, phi_c, phi_d)
!!       Implicit None
!!       Integer, Pointer:: theta_a(:), theta_b(:), theta_c(:)
!!       Integer, Pointer:: phi_a(:), phi_b(:), phi_c(:), phi_d(:)
!!     End Subroutine topology_generate_dihe
!!
!!   FUNCTION
!!     Using a list of dihe, generate a list of dihe
!!
!!   AUTHOR
!!     IKUO 08.01.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
!Use information from bond list to generate a list of dihedral angles
SUBROUTINE topology_generate_dihe(topology)
  IMPLICIT NONE

!Arguments
  TYPE ( topology_parameters_type ), INTENT (INOUT)   :: topology

!Locals
  INTEGER :: nphi
  INTEGER :: first,last,search
  INTEGER :: firstb,lastb,searchb
  INTEGER :: ntheta
  INTEGER :: itheta,iitheta,iiitheta,iphi
  INTEGER, POINTER :: tmp_a(:),tmp_b(:),tmp_c(:),tmp_d(:)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
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
!!     Subroutine topology_connectivity_pack(bond_a, bond_b, theta_a, theta_b, theta_c,&
!!         phi_a, phi_b, phi_c, phi_d, label_atmname, label_resname,&
!!         atm_mass, atm_charge, map_mol_typ, map_mol_num, molecule_set, setup,&
!!         globenv)
!!       Implicit None
!!       Type(molecule_type), Dimension(:), Pointer:: molecule_set
!!       Type(setup_parameters_type), Intent (IN):: setup
!!       Type(global_environment_type), Intent (IN):: globenv
!!       Character(Len=20), Pointer:: label_resname(:)
!!       Character(Len=20), Pointer:: label_atmname(:)
!!       Real(Kind=dbl), Pointer:: atm_charge(:)
!!       Real(Kind=dbl), Pointer:: atm_mass(:)
!!       Integer, Pointer:: map_mol_typ(:)
!!       Integer, Pointer:: map_mol_num(:)
!!       Integer, Pointer:: bond_a(:), bond_b(:)
!!       Integer, Pointer:: theta_a(:), theta_b(:), theta_c(:)
!!       Integer, Pointer:: phi_a(:), phi_b(:), phi_c(:), phi_d(:)
!!     End Subroutine topology_connectivity_pack
!!
!!   FUNCTION
!!     topology connectivity pack
!!
!!   AUTHOR
!!     IKUO 08.17.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE topology_connectivity_pack(molecule_kind_set,molecule_set,&
                                      topology,globenv)

  IMPLICIT NONE

! Arguments
  TYPE ( molecule_type ), DIMENSION (:), POINTER :: molecule_set
  TYPE ( molecule_kind_type), DIMENSION(:), POINTER :: molecule_kind_set
  TYPE ( topology_parameters_type ), INTENT (INOUT) :: topology
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Locals
  CHARACTER ( LEN=20 ) :: label,str2
  INTEGER :: iw,ierror,ilen,n
  INTEGER :: iatom,natom
  CHARACTER (LEN=20), POINTER :: atm_mol_name(:)
  INTEGER, POINTER :: atm_res_num(:)
  INTEGER :: ibond,nbond
  INTEGER :: itheta,ntheta
  INTEGER :: iphi,nphi
  INTEGER :: iunit,counter
  INTEGER :: i,j,k,first,last,search
  INTEGER :: nmol_typ,nmol_tot,nmol
  INTEGER, POINTER :: checkme(:)
  INTEGER, POINTER :: first_list(:)
  CHARACTER ( LEN=20 ) :: nhcopt
  TYPE (molecule_kind_type), POINTER :: molecule_kind
  TYPE (molecule_type), POINTER :: molecule
  TYPE(bond_type), DIMENSION(:), POINTER            :: bond_list
  TYPE(bend_type), DIMENSION(:), POINTER            :: bend_list
  TYPE(torsion_type), DIMENSION(:), POINTER            :: torsion_list
  INTEGER, DIMENSION(:), POINTER         :: molecule_list
  TYPE(atom_type), DIMENSION(:), POINTER            :: atom_list
  TYPE(local_molecule_type), POINTER                :: lmi

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  iw = globenv%scr
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  natom  = topology%natoms
  nmol_typ = topology%map_mol_typ(natom)
  IF(nmol_typ < 0) THEN
    CALL stop_program ( "topology_connectivity_pack",&
                        "no moleculular type defined")
  ELSE
    !Allocate space for the molecules_kind_set
    NULLIFY(molecule_kind_set)
    CALL allocate_molecule_kind_set(molecule_kind_set,nmol_typ)
  END IF
!------------------------------------------------------------------------------
  ALLOCATE(first_list(100000))
  nmol_tot=0
  counter=0
  DO i=1,nmol_typ
    ALLOCATE(molecule_list(100000))
    CALL find_boundary(topology%map_mol_typ,natom,first,last,i)
    nmol_tot=nmol_tot+topology%map_mol_num(last)
    DO j=1,topology%map_mol_num(last)
      counter=counter+1
      CALL find_boundary2(topology%map_mol_typ,topology%map_mol_num,natom,first,last,i,j)
      first_list(counter)=first
      molecule_list(counter)=counter
    END DO
    CALL reallocate(molecule_list,1,counter)
    molecule_kind => molecule_kind_set(i)
    CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           molecule_list=molecule_list)
    DEALLOCATE(molecule_list)
  END DO

  CALL allocate_molecule_set(molecule_set,nmol_tot)
  DO i=1,nmol_tot
    molecule => molecule_set(i)
    ALLOCATE(lmi,STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('topology_connectivity_pack','atom_list',natom)
    lmi%nstates=0
    NULLIFY(lmi%states)
    CALL set_molecule(molecule,lmi=lmi)
  END DO

  CALL reallocate(first_list,1,nmol_tot)
  CALL set_molecule_set(molecule_set,first_list)
  DEALLOCATE(first_list)

  DO i=1,nmol_typ
    natom  = SIZE(topology%map_mol_typ)
    nbond  = SIZE(topology%bond_a)
    ntheta = SIZE(topology%theta_a)
    nphi   = SIZE(topology%phi_a)
    nhcopt = "GLOBAL"
    molecule_kind => molecule_kind_set(i)

    CALL find_boundary(topology%map_mol_typ,natom,first,last,i)
    nmol=topology%map_mol_num(last)
    CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           kind_number=i,&
                           name=topology%label_resname(first),&
                           nhcopt=nhcopt)
    DO j=1,nmol
      molecule => molecule_set(j)
      CALL set_molecule(molecule,molecule_kind=molecule_kind)
      CALL get_molecule(molecule,molecule_kind)
    END DO

    natom = (last-first+1)/ topology%map_mol_num(last)
    ALLOCATE(atom_list(natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('topology_connectivity_pack','atom_list',natom)
    DO j=1,natom
      atom_list(j)%name=topology%label_atmname(j+first-1)
    END DO

    CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           atom_list=atom_list)

    !--------------------------------------------------------------------------
    !Need to be generalized!!!
    !Figure out the number of atoms and atom types-----------------------------
    !j = last - first + 1
    !j = j / map_mol_num(last)
    !molecule_set(i)%molpar%natom=j
    !molecule_set(i)%molpar%natom_type=1
    !DO j=1,molecule_set(i)%molpar%natom-1
    !  IF(.NOT.(label_atmname(first+j-1)==label_atmname(first+j))) THEN
    !    molecule_set(i)%molpar%natom_type=molecule_set(i)%molpar%natom_type + 1
    !  END IF
    !END DO
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !get the memory first
    !ALLOCATE(molecule_set(i)%molpar%aname(molecule_set(i)%molpar%natom),STAT=ierror)
    !IF(ierror/=0) CALL stop_memory &
    ! ('topology_connectivity_pack','aname',molecule_set(i)%molpar%natom)
    !LLOCATE(molecule_set(i)%molpar%alabel(molecule_set(i)%molpar%natom),STAT=ierror)
    !F(ierror/=0) CALL stop_memory &
    ! ('topology_connectivity_pack','alabel',molecule_set(i)%molpar%natom)
    !LLOCATE(molecule_set(i)%molpar%aweight(molecule_set(i)%molpar%natom),STAT=ierror)
    !F(ierror/=0) CALL stop_memory &
    ! ('topology_connectivity_pack','aweight',molecule_set(i)%molpar%natom)
    !LLOCATE(molecule_set(i)%molpar%acharge(molecule_set(i)%molpar%natom),STAT=ierror)
    !F(ierror/=0) CALL stop_memory &
    ! ('topology_connectivity_pack','acharge',molecule_set(i)%molpar%natom)
    !O j=1,molecule_set(i)%molpar%natom
    ! molecule_set(i)%molpar%aname(j)   = ''
    ! molecule_set(i)%molpar%aweight(j) = 0.0
    ! molecule_set(i)%molpar%acharge(j) = 0.0
      !Stuff the aname,alabel,aweight,acharge
    ! molecule_set(i)%molpar%aname(j)   = label_atmname(first+j-1)
    ! molecule_set(i)%molpar%aweight(j) = atm_mass(first+j-1)
    ! CALL convert_to_cp2k_units(mass = molecule_set(i)%molpar%aweight(j))
    ! molecule_set(i)%molpar%acharge(j) = atm_charge(first+j-1)
    !ND DO
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Stuff the bond info
    ALLOCATE(checkme(nbond),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('topology_connectivity_pack','checkme',nbond)
    checkme(:)=0
    DO j=1,natom
      IF ((topology%map_mol_num(j)==1).AND.(topology%map_mol_typ(j)==i)) THEN
        DO k=1,nbond
          IF ((topology%bond_a(k)==j).OR.(topology%bond_b(k)==j)) THEN
            checkme(k)=1
          END IF
        END DO
      END IF
    END DO
    nbond = SUM(checkme)
    ALLOCATE(bond_list(nbond),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('topology_connectivity_pack','bond_list',nbond)
    k=0
    DO j=1,nbond
      IF (checkme(j)==1) THEN
        k=k+1
        !-------------------------------------------------------------------
        !Check and make sure
        bond_list(i)%a=topology%bond_a(j)
        bond_list(i)%b=topology%bond_b(j)
        !molecule_set(i)%molpar%bonds(3,k) = j
        !point this to the right bond_kind_type
        !-------------------------------------------------------------------
      END IF
    END DO
    IF (k /= nbond) THEN
      !ERROR STOP PROGRAM
      CALL stop_program ( "topology_connectivity_pack", "nbond counting error")
    END IF
    DEALLOCATE(checkme,STAT=ierror)
    IF (ierror /= 0) CALL stop_memory ('read_molecule_section','checkme')
    CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           bond_list=bond_list,&
                           nbond=nbond)



    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Stuff the theata info
    ALLOCATE(checkme(ntheta),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('topology_connectivity_pack','checkme',ntheta)
    checkme(:)=0
    DO j=1,natom
      IF ((topology%map_mol_num(j)==1).AND.(topology%map_mol_typ(j)==i)) THEN
        DO k=1,ntheta
          IF ((topology%theta_a(k)==j).OR.(topology%theta_b(k)==j).OR.(topology%theta_c(k)==j)) THEN
            checkme(k)=1
          END IF
        END DO
      END IF
    END DO
    ntheta = SUM(checkme)
    ALLOCATE(bend_list(ntheta),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('topology_connectivity_pack','bend_list',ntheta)
    k=0
    DO j=1,ntheta
      IF (checkme(j)==1) THEN
        k=k+1
        !-------------------------------------------------------------------
        !Check and make sure
        bend_list(i)%a=topology%theta_a(j)
        bend_list(i)%b=topology%theta_b(j)
        bend_list(i)%c=topology%theta_c(j)
        !molecule_set(i)%molpar%bends(4,k) = j
        !point this to the right bend_kind_type
        !-------------------------------------------------------------------
      END IF
    END DO
    IF (k /= ntheta) THEN
      !ERROR STOP PROGRAM
      CALL stop_program ( "topology_connectivity_pack", "nbends counting error")
    END IF
    DEALLOCATE(checkme,STAT=ierror)
    IF (ierror /= 0) CALL stop_memory ('read_molecule_section','checkme')
    CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           bend_list=bend_list,&
                           nbend=ntheta)



    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Stuff the torsion info
    ALLOCATE(checkme(nphi),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('topology_connectivity_pack','checkme',nphi)
    checkme(:)=0
    DO j=1,natom
      IF ((topology%map_mol_num(j)==1).AND.(topology%map_mol_typ(j)==i)) THEN
        DO k=1,nphi
          IF ((topology%phi_a(k)==j).OR.(topology%phi_b(k)==j).OR.(topology%phi_c(k)==j)&
              .OR.(topology%phi_d(k)==j)) THEN
            checkme(k)=1
          END IF
        END DO
      END IF
    END DO
    nphi = SUM(checkme)
    ALLOCATE(torsion_list(nphi),STAT=ierror)
    k=0
    DO j=1,nphi
      IF (checkme(j)==1) THEN
        k=k+1
        !-------------------------------------------------------------------
        !Check and make sure
        torsion_list(i)%a=topology%phi_a(j)
        torsion_list(i)%b=topology%phi_b(j)
        torsion_list(i)%c=topology%phi_c(j)
        torsion_list(i)%d=topology%phi_d(j)
        !molecule_set(i)%molpar%torsions(5,k) = j
        !point this to the right torsion_kind_type
        !-------------------------------------------------------------------
      END IF
    END DO
    IF (k /= nphi) THEN
      !ERROR STOP PROGRAM
      CALL stop_program ( "topology_connectivity_pack", "ntorsions counting error")
    END IF
    DEALLOCATE(checkme,STAT=ierror)
    IF (ierror /= 0) CALL stop_memory ('read_molecule_section','checkme')
    CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           torsion_list=torsion_list,&
                           ntorsion=nphi)
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
  END DO

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
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
!!     Subroutine topology_coordinate_pack(r, label_atmname, particle_set,&
!!         atomic_kind_set, molecule_set, setup, globenv)
!!       Implicit None
!!       Real(Kind=dbl), Dimension(:,:), Pointer:: r
!!       Character(Len=20), Pointer:: label_atmname(:)
!!       Type(global_environment_type), Intent (IN):: globenv
!!       Type(atomic_kind_type), Dimension(:), Pointer:: atomic_kind_set
!!       Type(particle_type), Dimension(:), Pointer:: particle_set
!!       Type(molecule_type), Dimension(:), Intent (IN):: molecule_set
!!       Type(setup_parameters_type):: setup
!!     End Subroutine topology_coordinate_pack
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

  IMPLICIT NONE

! Arguments
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
  TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
  TYPE ( molecule_type ), DIMENSION ( : ) :: molecule_set
  TYPE (molecule_kind_type), DIMENSION (:),POINTER::molecule_kind_set
  TYPE ( topology_parameters_type ), INTENT (INOUT) :: topology

! Locals
  INTEGER :: ierror, ilen, iw,  ia, ie, i,nmol_typ
  INTEGER, PARAMETER :: max_name_length = 60
  INTEGER, PARAMETER :: nblock_atom = 1000,&
                        nblock_kind = 10

  TYPE(atomic_kind_type), POINTER :: atomic_kind
  TYPE(epc_potential_type), POINTER :: epc_potential
  CHARACTER(LEN=max_name_length)  :: current_label,label
  CHARACTER(LEN=200)              :: message
  CHARACTER(LEN=8)                :: unit_name
  CHARACTER(LEN=3)                :: test_result
  CHARACTER(LEN=20),DIMENSION(:),ALLOCATABLE       :: tuple,work
  CHARACTER ( LEN = 80 ) :: bl = ' '

  INTEGER                         :: iatom,ikind,istat,natom,nkind,z,n,iat,&
                                     natom_types,nmol_type,k,j
  LOGICAL                         :: new_kind,coords,&
                                       init
  INTEGER :: imol, na, nat, j1, j2, nf, nt, i1, i2, i3, i4
  INTEGER :: isos, allgrp, npol, nrho0, anum,ios
  CHARACTER ( LEN = 20 ) :: libstring, string, string2, at1, at2, at3, at4


  INTEGER, DIMENSION(:), POINTER :: kind_of,natom_of_kind,natom_of_kind_check
  INTEGER, DIMENSION(:,:), POINTER  :: atomlist,kind_label
  REAL(dbl), DIMENSION (:), POINTER :: charge, mass

  TYPE (molecule_kind_type), POINTER ::molecule_kind
  INTEGER :: nmolecule
  TYPE(atom_type), DIMENSION(:), POINTER            :: atom_list


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  iw = globenv % scr

  nkind = 0
  natom = SIZE(topology%label_atmname)
  iat = 0

! count total number of molecules and atoms
  nmol_type = SIZE(molecule_kind_set)
  topology % nmol = 0
  topology % natoms = 0
  topology % natom_type = 0
  n = 0
  ALLOCATE ( work( 10000 ), STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'read_molecule_section', &
       'work', topology % natom_type )
  DO i = 1, nmol_type
     molecule_kind => molecule_kind_set(i)
     CALL get_molecule_kind(molecule_kind=molecule_kind,&
                            nmolecule=nmolecule,&
                            natom=natom,&
                            atom_list=atom_list)

     topology%nmol = topology % nmol + nmolecule
     topology%natoms = topology % natoms + nmolecule * natom
     !topology%natom_type = topology % natom_type +  &
     !     molecule_set(i) % molpar % natom_type
     DO j = 1, natom

        IF(str_search(work,n,atom_list(j)%name)==0)THEN
           n = n + 1
           work(n) = atom_list(j)%name
        END IF
     END DO
  END DO
  DEALLOCATE(work)
  topology%natom_type = n

  ALLOCATE(mass(topology%natom_type),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','mass')
! generate unique atom list
  ALLOCATE ( work( topology % natom_type ), STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'read_molecule_section', &
       'work', topology % natom_type )

  n = 0
  nt = 0
  DO i = 1, nmol_type
     molecule_kind => molecule_kind_set(i)
     CALL get_molecule_kind(molecule_kind=molecule_kind,&
                            nmolecule=nmolecule,&
                            natom=natom,&
                            atom_list=atom_list)
     DO j = 1, natom
        IF(str_search(work,n,atom_list(j)%name)==0)THEN
           n = n + 1
           work(n) = atom_list(j)%name
           mass(n) = topology%atm_mass(nt+j)
        END IF
     END DO
     nt = nt + natom
  END DO

  ALLOCATE ( topology%atom_names(n), STAT=ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'read_molecule_section', &
       'setup%atom_names', n )
  topology % natom_type = n
  topology % atom_names(1:n) = work(1:n)
  DEALLOCATE( work, STAT=ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'read_molecule_section', &
       'work' )
  
  CALL write_molecule_kind_set(molecule_kind_set,globenv)
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!Put all the info from temporary storage info compatible FIST storage

  natom = SIZE(topology%label_atmname)
  ALLOCATE(kind_label(max_name_length,natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','kind_label')
  ALLOCATE(kind_of(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','kind_of')
  ALLOCATE(natom_of_kind(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','natom_of_kind')

  kind_of(:) = 0
  natom_of_kind(:) = 0

  !Determine the number of kinds (nkind)
  nkind = 0
  DO iatom=1,natom
    new_kind = .TRUE.
    DO ikind=1,nkind
      CALL integer_to_string(kind_label(:,ikind),string)
      IF(string(1:20)==topology%label_atmname(iatom)) THEN
        new_kind = .FALSE.
        kind_of(iatom)=ikind
        natom_of_kind(ikind) = natom_of_kind(ikind)+1
        EXIT
      END IF
    END DO
    IF(new_kind) THEN
      nkind = nkind + 1
      kind_of(iatom)=nkind
      CALL string_to_integer(topology%label_atmname(iatom),kind_label(:,nkind))
      natom_of_kind(ikind) = natom_of_kind(ikind)+1
    END IF
  END DO
  CALL reallocate(kind_label,1,max_name_length,1,nkind)
  CALL reallocate(natom_of_kind,1,nkind)
  ALLOCATE(charge(nkind),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','charge')

  !Stuff the atomlist(1:nkind,1:natom_of_kind(ikind))
  NULLIFY (atomlist)
  ALLOCATE(atomlist(nkind,MAXVAL(natom_of_kind)),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','atomlist')
  DO ikind=1,nkind
    iat = 0
    DO iatom = 1,natom
      IF (kind_of(iatom)==ikind) THEN
        iat = iat + 1
        atomlist(ikind,iat) = iatom
      END IF
    END DO
  END DO


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

!   *** Allocate the data structure for the atomic kind information ***

  NULLIFY(atomic_kind_set)
  CALL allocate_atomic_kind_set(atomic_kind_set,nkind)

!   *** Initialize the data structure for atomic kind information ***

  DO ikind=1,nkind
    atomic_kind => atomic_kind_set(ikind)
    CALL integer_to_string(kind_label(1:20,ikind),string)
    CALL set_atomic_kind(atomic_kind=atomic_kind,&
                         kind_number=ikind,&
                         name=string,&
                         atom_list=atomlist(ikind,1:natom_of_kind(ikind)))
  END DO


  DO i=1,nmol_type
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           natom=natom,&
                           atom_list=atom_list)
    DO j=1,natom
      DO ikind=1,nkind
      !DO ikind=1,2
        atomic_kind =>atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                          name=label) 
        IF(atom_list(j)%name == label) THEN
          atom_list(j)%atomic_kind => atomic_kind_set(ikind)
        END IF
      END DO
    END DO
  END DO


!   *** Allocate the data structure for the atomic information ***
  natom = SIZE(topology%label_atmname)
  NULLIFY(particle_set)
  CALL allocate_particle_set(particle_set,natom)

!   *** Initialize the data structure for the atomic information ***
  DO iatom=1,natom
    ikind = kind_of(iatom)
    particle_set(iatom) % abase = 0
    particle_set(iatom) % grid_id = 0
    particle_set(iatom)%atomic_kind => atomic_kind_set(ikind)
    particle_set(iatom)%iatom = iatom
    particle_set(iatom)%f(:) = 0.0_dbl
    particle_set(iatom)%r(:) = topology%r(:,iatom)
    particle_set(iatom)%v(:) = 0.0_dbl
    NULLIFY ( particle_set ( iatom ) % coef_list )
  END DO

!   *** check if the labels in &COORD match with the &MOLECULE setup
!   *** and store the charge and mass read in the &MOLECULE section

  !CALL stop_program("topology_coordinate_pack","checkpoint")
!  ALLOCATE (natom_of_kind_check(nkind), STAT=istat) 
!  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','natom_of_kind_check')
!  natom_of_kind_check = 0
!  nmol_type = SIZE (molecule_set) 
!  natom_types = SIZE(atomic_kind_set)
!  natom_types = SIZE(topology%atom_names)
!  iat = 0
!  DO i = 1, nmol_type
!    DO k = 1, molecule_kind_set(i) % nmolecule
!      DO j = 1, molecule_kind_set(i) % natom
!        iat = iat + 1
!        n = str_search ( topology%atom_names, natom_types, &
!             molecule_set ( i ) % molpar % aname ( j ) )
!        charge(n)= molecule_set(i) % molpar % acharge(j)
!!dbg this should be done using general unit conversion routine
!!AU
!!            mass(n)= massunit* molecule_set(i) % molpar % aweight(j)
!! kelvin
!        mass(n)= molecule_set(i) % molpar % aweight(j)
!! gt
!        CALL get_atomic_kind(particle_set(iat)%atomic_kind,&
!                             kind_number=ikind)
!        IF (ikind/=n) CALL stop_program ( 'topology_coordinate_pack', &
!             'kind label ordering in coord and mol dont match ')
!        CALL integer_to_string(kind_label(1:20,ikind),string)
!        IF (string(1:20)/=molecule_set(i)%molpar%aname(j)) CALL stop_program ( &
!             'topology_coordinate_pack', &
!             'inconsistent kind labels in coord and mol')
!        natom_of_kind_check(ikind) = natom_of_kind_check(ikind) + 1
!        iatom =  atomlist(ikind,natom_of_kind_check(ikind)) 
!        IF (iatom/=iat) CALL stop_program ( 'topology_coordinate_pack', &
!             'atom lists in coord and mol dont match ')
!      END DO
!    END DO
!  END DO

  DO ikind=1,nkind
    atomic_kind => atomic_kind_set(ikind)
    CALL set_atomic_kind(atomic_kind=atomic_kind,&
                         mass=mass(ikind)) 
    CALL write_atomic_kind(atomic_kind,globenv)
  END DO

  IF (globenv%program_name == "FIST") THEN
    CALL  init_atomic_kind_set(atomic_kind_set,globenv)
    DO ikind=1,nkind
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                        epc_potential=epc_potential) 
      CALL set_potential(potential=epc_potential,&
                        qeff=charge(ikind))
    END DO
  END IF

!   *** Release work storage ***
  DEALLOCATE (atomlist,kind_label,kind_of,natom_of_kind,STAT=istat)
  IF (istat /= 0) THEN
    CALL stop_memory('topology_coordinate_pack',&
              'release work storage')
  END IF
  DEALLOCATE (charge,mass,STAT=istat)
  IF (istat /= 0) THEN
    CALL stop_memory('topology_coordinate_pack',&
              'release work storage')
  END IF

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  !CALL stop_program("topology_coordinate_pack","checkpoint end")

END SUBROUTINE topology_coordinate_pack


END MODULE topology_util
