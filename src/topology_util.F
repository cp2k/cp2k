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
                                             set_atomic_kind
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE external_potential_types,        ONLY: epc_potential_type,&
                                             set_potential
  USE global_types,                    ONLY: global_environment_type
  USE input_types,                     ONLY: setup_parameters_type
  USE kinds,                           ONLY: dbl
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_bcast
  USE molecule_types,                  ONLY: init_molecule_type,&
                                             molecule_type
  USE nrutil,                          ONLY: swap
  USE particle_types,                  ONLY: allocate_particle_set,&
                                             particle_type
  USE periodic_table,                  ONLY: ptable
  USE simulation_cell,                 ONLY: cell_type,&
                                             pbc
  USE string_utilities,                ONLY: integer_to_string,&
                                             str_search,&
                                             string_to_integer,&
                                             uppercase,&
                                             xstring
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  PRIVATE
  PUBLIC :: topology_generate_molecule,azero,find_boundary,topology_coord_pack,&
            topology_generate_bond,topology_generate_bend,&
            topology_generate_dihe,sort_bond,&
            sort_bend,sort_dihe,topology_pack

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
SUBROUTINE topology_generate_molecule(natom,atm_mol_name,nbond,bond_a,bond_b,&
                                      map_mol_typ,map_mol_num)

  IMPLICIT NONE

! Arguments
  INTEGER :: natom,nbond
  CHARACTER (LEN=20), POINTER :: atm_mol_name(:)
  INTEGER, POINTER :: bond_a(:)
  INTEGER, POINTER :: bond_b(:)
  INTEGER, POINTER :: map_mol_typ(:) !map_mol_typ(iatom)= mol_typ number
  INTEGER, POINTER :: map_mol_num(:) !map_mol_num(iatom)= mol number

! Locals
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

  ALLOCATE(checkme(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_generate_molecule','checkme')
  !Zero the arrays
  CALL azero(checkme,natom)
  CALL azero(map_mol_typ,natom)
  CALL azero(map_mol_num,natom)

  !Parse the atoms list to see how many different molecule types there are
  ntype = 1
  map_mol_typ(1)=1
  DO iatom=2,natom
    IF(atm_mol_name(iatom-1)==atm_mol_name(iatom)) THEN
      map_mol_typ(iatom)=ntype
    ELSE
      ntype=ntype+1
      map_mol_typ(iatom)=ntype
    END IF
  END DO

  DO itype=1,ntype
    !search for first atom in the type list and last to cut down time
    CALL find_boundary(map_mol_typ,natom,first,last,itype)

    inum = 0

    DO iatom=first,last
      CALL azero(checkme,natom)

      IF(map_mol_num(iatom)==0) THEN

        inum=inum+1

        checkme(iatom)=1
       
        DO WHILE (SUM(checkme)>0)
          DO iiatom=first,last
            IF((checkme(iiatom)==1)) THEN
              checkme(iiatom)=0
              map_mol_num(iiatom)=inum
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
                IF (bond_a(iii)==iiatom) THEN
                  checkme(bond_b(iii))=1
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
!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     azero
!!
!!   SYNOPSIS
!!     Subroutine azero(num_array, ntot)
!!       Implicit None
!!       Integer, Pointer:: num_array(:)
!!       Integer:: ntot
!!     End Subroutine azero
!!
!!   FUNCTION
!!     zero the sucker
!!
!!   AUTHOR
!!     IKUO 08.01.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE azero(num_array,ntot)

  IMPLICIT NONE

! Arguments
  INTEGER, POINTER :: num_array(:)
  INTEGER :: ntot

! Locals
  INTEGER :: i

!------------------------------------------------------------------------------
  DO i=1,ntot
    num_array(i) = 0
  END DO

END SUBROUTINE azero


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
SUBROUTINE topology_generate_bond(cell,r,label_atmname,bond_a,bond_b,&
                                  atm_mass,atm_charge)
  IMPLICIT NONE
!Arguments
  TYPE ( cell_type ), POINTER :: cell
  REAL(dbl), DIMENSION(:,:), POINTER :: r
  CHARACTER (LEN=20), POINTER :: label_atmname(:)
  INTEGER, POINTER :: bond_a(:),bond_b(:)
  REAL (dbl), POINTER :: atm_charge(:)
  REAL (dbl), POINTER :: atm_mass(:)

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
  natom = SIZE(r(1,:))

  !Get the temporary array that we need
  ALLOCATE(tbond_a(2*natom))
  ALLOCATE(tbond_b(2*natom))
  ALLOCATE(pos_x(natom))
  ALLOCATE(pos_y(natom))
  ALLOCATE(pos_z(natom))
  ALLOCATE(atm_symbol(natom))
  ALLOCATE(covalent_radius(natom))
  ALLOCATE(atm_mass(natom))
  ALLOCATE(atm_charge(natom))
!------------------------------------------------------------------------------
  DO iatom=1,natom
    atm_symbol(iatom) = label_atmname(iatom)
  END DO
!------------------------------------------------------------------------------
  WRITE(*,*) "GENERATING BONDS"
  !Set up the arrays that we need for this routine
  CALL azero(tbond_a,2*natom)
  CALL azero(tbond_b,2*natom)

  DO iatom=1,natom
    pos_x(iatom) = r(1,iatom)
    pos_y(iatom) = r(2,iatom)
    pos_z(iatom) = r(3,iatom)
    
    nelem = SIZE(ptable)
    DO ielem=1,nelem 
      IF(ptable(ielem)%symbol==atm_symbol(iatom)) THEN
        covalent_radius(iatom) = ptable(ielem)%covalent_radius
        CALL convert_to_cp2k_units("ANGSTROM",length=covalent_radius(iatom))
        atm_mass(iatom) = ptable(ielem)% amass
        !Need to figure out the correct charge!!!
        atm_charge(iatom) = 0.0
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
      rr = pbc ( dr, cell )
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
      END IF
    END DO
  END DO
  WRITE(*,*) "NBOND = ",nbond
!------------------------------------------------------------------------------
  !Copy the stuff from temporary array into the more permanent one
  ALLOCATE(bond_a(nbond))
  ALLOCATE(bond_b(nbond))
  DO ibond=1,nbond
    bond_a(ibond) = tbond_a(ibond)
    bond_b(ibond) = tbond_b(ibond)
  END DO
!------------------------------------------------------------------------------
  !Sort out the bond list
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
!------------------------------------------------------------------------------
  !Sort out the duplicates... (Shouldn't have any but just to check)
  DO ibond=2,nbond
    IF((bond_a(ibond-1)==bond_a(ibond)).AND.&
       (bond_b(ibond-1)==bond_b(ibond))) THEN
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
SUBROUTINE topology_generate_bend(bond_a,bond_b,&
                                   theta_a,theta_b,theta_c)
  IMPLICIT NONE

!Arguments
  INTEGER, POINTER :: bond_a(:),bond_b(:)
  INTEGER :: ntheta
  INTEGER, POINTER :: theta_a(:),theta_b(:),theta_c(:)

!Locals
  INTEGER :: first,last,search,countme
  INTEGER :: nbond
  INTEGER :: ibond,iibond,itheta
  INTEGER, POINTER :: tmp_a(:),tmp_b(:),tmp_c(:)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  nbond = SIZE(bond_a)
  ntheta = 0
  !Stupid and inefficient but faster than alternative
  ALLOCATE(tmp_a(2*nbond))
  ALLOCATE(tmp_b(2*nbond))
  ALLOCATE(tmp_c(2*nbond))
!------------------------------------------------------------------------------
  ntheta=0
  DO ibond=1,nbond
    !Case A
    search = bond_a(ibond)
    CALL find_boundary(bond_a,nbond,first,last,search)
    IF (first /= 0) THEN
      DO iibond=first,last
        IF (bond_b(ibond)/=bond_b(iibond)) THEN
          ntheta = ntheta + 1
          tmp_a(ntheta) = bond_b(ibond)
          tmp_b(ntheta) = search
          tmp_c(ntheta) = bond_b(iibond)
        END IF
      END DO
    END IF
    !Case B
    search = bond_b(ibond)
    CALL find_boundary(bond_a,nbond,first,last,search)
    IF (first /= 0) THEN
      DO iibond=first,last
        IF(bond_a(ibond)/=bond_b(iibond)) THEN
          ntheta = ntheta + 1
          tmp_a(ntheta) = bond_a(ibond)
          tmp_b(ntheta) = search
          tmp_c(ntheta) = bond_b(iibond)
        END IF
      END DO
    END IF
    !Case C
    search = bond_b(ibond)
    CALL find_boundary(bond_b,nbond,first,last,search)
    IF (first /= 0) THEN
      DO iibond=first,last
        IF (bond_a(ibond)/=bond_a(iibond)) THEN
          ntheta = ntheta + 1
          tmp_a(ntheta) = bond_a(ibond)
          tmp_b(ntheta) = search
          tmp_c(ntheta) = bond_a(iibond)
        END IF
      END DO
    END IF
  END DO
  WRITE(*,*) "NTHETA = ",ntheta
!------------------------------------------------------------------------------
  !Sort all the bends first
  CALL sort_bend(ntheta,tmp_a,tmp_b,tmp_c)

  ALLOCATE(theta_a(ntheta))
  ALLOCATE(theta_b(ntheta))
  ALLOCATE(theta_c(ntheta))
  countme=1
  theta_a(countme) = tmp_a(1)
  theta_b(countme) = tmp_b(1)
  theta_c(countme) = tmp_c(1)
  DO itheta=2,ntheta
    IF( .NOT.(tmp_a(itheta) == theta_a(countme) .AND. &
        tmp_b(itheta) == theta_b(countme) .AND. &
        tmp_c(itheta) == theta_c(countme))) THEN
      countme = countme + 1
      theta_a(countme) = tmp_a(itheta)
      theta_b(countme) = tmp_b(itheta)
      theta_c(countme) = tmp_c(itheta)
    END IF
  END DO
  ntheta = countme
  WRITE(*,*) "NTHETA = ",ntheta

  CALL REALLOCATE(theta_a,1,ntheta)
  CALL REALLOCATE(theta_b,1,ntheta)
  CALL REALLOCATE(theta_c,1,ntheta)
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
SUBROUTINE topology_generate_dihe(theta_a,theta_b,theta_c,&
                                  phi_a,phi_b,phi_c,phi_d)
  IMPLICIT NONE

!Arguments
  INTEGER, POINTER :: theta_a(:),theta_b(:),theta_c(:)
  INTEGER, POINTER :: phi_a(:),phi_b(:),phi_c(:),phi_d(:)

!Locals
  INTEGER :: nphi
  INTEGER :: first,last,search
  INTEGER :: firstb,lastb,searchb
  INTEGER :: ntheta
  INTEGER :: itheta,iitheta,iiitheta,iphi
  INTEGER, POINTER :: tmp_a(:),tmp_b(:),tmp_c(:),tmp_d(:)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  ntheta = SIZE(theta_a)
  nphi = 0
  !Stupid and inefficient but faster than alternative
  ALLOCATE(tmp_a(2*ntheta))
  ALLOCATE(tmp_b(2*ntheta))
  ALLOCATE(tmp_c(2*ntheta))
  ALLOCATE(tmp_d(2*ntheta))
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  ALLOCATE(phi_a(nphi))
  ALLOCATE(phi_b(nphi))
  ALLOCATE(phi_c(nphi))
  ALLOCATE(phi_d(nphi))
  DO iphi=1,nphi
    phi_a(iphi) = tmp_a(iphi)
    phi_b(iphi) = tmp_b(iphi)
    phi_c(iphi) = tmp_c(iphi)
    phi_d(iphi) = tmp_d(iphi)
  END DO
  nphi = SIZE(phi_a)
  WRITE(*,*) "NPHI = ",nphi
!------------------------------------------------------------------------------
  DEALLOCATE(tmp_a)
  DEALLOCATE(tmp_b)
  DEALLOCATE(tmp_c)
  DEALLOCATE(tmp_d)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE topology_generate_dihe
 


!******************************************************************************
!!****** topology_util/topology_coord_pack [1.0] *
!!
!!   NAME
!!     topology_coord_pack
!!
!!   SYNOPSIS
!!     Subroutine topology_coord_pack(r, label_atmname, particle_set,&
!!         atomic_kind_set, mol_setup, setup, globenv)
!!       Implicit None
!!       Real(Kind=dbl), Dimension(:,:), Pointer:: r
!!       Character(Len=20), Pointer:: label_atmname(:)
!!       Type(global_environment_type), Intent (IN):: globenv
!!       Type(atomic_kind_type), Dimension(:), Pointer:: atomic_kind_set
!!       Type(particle_type), Dimension(:), Pointer:: particle_set
!!       Type(molecule_type), Dimension(:), Intent (IN):: mol_setup
!!       Type(setup_parameters_type):: setup
!!     End Subroutine topology_coord_pack
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

SUBROUTINE topology_coord_pack (r,label_atmname,particle_set, &
                              atomic_kind_set, &
                              mol_setup,  setup, globenv )

  IMPLICIT NONE

! Arguments
  REAL(dbl), DIMENSION(:,:), POINTER :: r
  CHARACTER (LEN=20), POINTER :: label_atmname(:)
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
  TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
  TYPE ( molecule_type ), DIMENSION ( : ), INTENT ( IN ) :: mol_setup
  TYPE ( setup_parameters_type ) :: setup

! Locals
  INTEGER :: ierror, ilen, iw,  ia, ie, i
  INTEGER, PARAMETER :: max_name_length = 60
  INTEGER, PARAMETER :: nblock_atom = 1000,&
                        nblock_kind = 10

  TYPE(atomic_kind_type), POINTER :: atomic_kind
  TYPE(epc_potential_type), POINTER :: epc_potential
  CHARACTER(LEN=max_name_length)  :: current_label,label,string
  CHARACTER(LEN=200)              :: message
  CHARACTER(LEN=8)                :: unit_name
  CHARACTER(LEN=3)                :: test_result
  INTEGER                         :: iatom,ikind,istat,natom,nkind,z,n,iat,&
                                     natom_types,nmol_type,k,j
  LOGICAL                         :: new_kind,coords,&
                                       init

  INTEGER, DIMENSION(:), POINTER :: kind_of,natom_of_kind,natom_of_kind_check
  INTEGER, DIMENSION(:,:), POINTER  :: atom_list,kind_label
  REAL(dbl), DIMENSION (:), POINTER :: charge, mass

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  iw = globenv % scr

  nkind = 0
  natom = SIZE(label_atmname)
  WRITE(iw,*) "NATOM = ",natom
  iat = 0
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!Put all the info from temporary storage info compatible FIST storage

  ALLOCATE(kind_label(max_name_length,natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coord_pack','kind_label')
  ALLOCATE(kind_of(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coord_pack','kind_of')
  ALLOCATE(natom_of_kind(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coord_pack','natom_of_kind')
  CALL azero(kind_of,natom)
  CALL azero(natom_of_kind,natom)

  !Determine the number of kinds (nkind)
  nkind = 0
  DO iatom=1,natom
    new_kind = .TRUE.
    DO ikind=1,nkind
      CALL integer_to_string(kind_label(:,ikind),string)
      IF(string(1:20)==label_atmname(iatom)) THEN
        new_kind = .FALSE.
        kind_of(iatom)=ikind
        natom_of_kind(ikind) = natom_of_kind(ikind)+1
        EXIT
      END IF
    END DO
    IF(new_kind) THEN
      nkind = nkind + 1
      kind_of(iatom)=nkind
      CALL string_to_integer(label_atmname(iatom),kind_label(:,nkind))
      natom_of_kind(ikind) = natom_of_kind(ikind)+1
    END IF
  END DO
  WRITE(iw,*) "NKIND = ",nkind
  CALL reallocate(kind_label,1,max_name_length,1,nkind)
  CALL reallocate(natom_of_kind,1,nkind)
  ALLOCATE(charge(nkind),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coord_pack','charge')
  ALLOCATE(mass(nkind),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coord_pack','charge')

  !Stuff the atom_list(1:nkind,1:natom_of_kind(ikind))
  NULLIFY (atom_list)
  ALLOCATE(atom_list(nkind,MAXVAL(natom_of_kind)),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coord_pack','atom_list')
  DO ikind=1,nkind
    iat = 0
    DO iatom = 1,natom
      IF (kind_of(iatom)==ikind) THEN
        iat = iat + 1
        atom_list(ikind,iat) = iatom
      END IF
    END DO
  END DO

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

!   *** Allocate the data structure for the atomic kind information ***

  CALL allocate_atomic_kind_set(atomic_kind_set,nkind)

!   *** Initialize the data structure for atomic kind information ***

  DO ikind=1,nkind
    atomic_kind => atomic_kind_set(ikind)
    CALL integer_to_string(kind_label(1:20,ikind),string)
    WRITE(*,*) "kind_label ",string(1:20)
    CALL set_atomic_kind(atomic_kind=atomic_kind,&
                         kind_number=ikind,&
                         name=string,&
                         atom_list=atom_list(ikind,1:natom_of_kind(ikind)))
  END DO

!   *** Allocate the data structure for the atomic information ***
  CALL allocate_particle_set(particle_set,natom)


!   *** Initialize the data structure for the atomic information ***
  DO iatom=1,natom
    ikind = kind_of(iatom)
    particle_set(iatom) % abase = 0
    particle_set(iatom) % grid_id = 0
    particle_set(iatom)%atomic_kind => atomic_kind_set(ikind)
    particle_set(iatom)%iatom = iatom
    particle_set(iatom)%f(:) = 0.0_dbl
    particle_set(iatom)%r(:) = r(:,iatom)
    particle_set(iatom)%v(:) = 0.0_dbl
    NULLIFY ( particle_set ( iatom ) % coef_list )
  END DO

!   *** check if the labels in &COORD match with the &MOLECULE setup
!   *** and store the charge and mass read in the &MOLECULE section

  !CALL stop_program("topology_coord_pack","checkpoint")
  ALLOCATE (natom_of_kind_check(nkind), STAT=istat) 
  IF (istat /= 0) CALL stop_memory('topology_coord_pack','natom_of_kind_check')
  natom_of_kind_check = 0
  nmol_type = SIZE (mol_setup) 
  natom_types = SIZE(setup%atom_names)
  WRITE(*,*) "nmol_type = ",nmol_type
  WRITE(*,*) "natom_type = ",natom_types
  iat = 0
  DO i = 1, nmol_type
    DO k = 1, mol_setup(i) % num_mol
      DO j = 1, mol_setup(i) % molpar % natom
        iat = iat + 1
        n = str_search ( setup%atom_names, natom_types, &
             mol_setup ( i ) % molpar % aname ( j ) )
        charge(n)= mol_setup(i) % molpar % acharge(j)
!dbg this should be done using general unit conversion routine
!AU
!            mass(n)= massunit* mol_setup(i) % molpar % aweight(j)
! kelvin
        mass(n)= mol_setup(i) % molpar % aweight(j)
! gt
        CALL get_atomic_kind(particle_set(iat)%atomic_kind,&
                             kind_number=ikind)
        IF (ikind/=n) CALL stop_program ( 'topology_coord_pack', &
             'kind label ordering in coord and mol dont match ')
        CALL integer_to_string(kind_label(1:20,ikind),string)
        IF (string(1:20)/=mol_setup(i)%molpar%aname(j)) CALL stop_program ( &
             'topology_coord_pack', &
             'inconsistent kind labels in coord and mol')
        natom_of_kind_check(ikind) = natom_of_kind_check(ikind) + 1
        iatom =  atom_list(ikind,natom_of_kind_check(ikind)) 
        IF (iatom/=iat) CALL stop_program ( 'topology_coord_pack', &
             'atom lists in coord and mol dont match ')
      END DO
    END DO
  END DO

  DO ikind=1,nkind
    atomic_kind => atomic_kind_set(ikind)
    CALL set_atomic_kind(atomic_kind=atomic_kind,&
                         mass=mass(ikind)) 
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
  DEALLOCATE (atom_list,kind_label,kind_of,natom_of_kind,&
              natom_of_kind_check,r,STAT=istat)
  IF (istat /= 0) THEN
    CALL stop_memory('topology_coord_pack',&
              'release work storage')
  END IF
  DEALLOCATE (charge,mass,STAT=istat)
  IF (istat /= 0) THEN
    CALL stop_memory('topology_coord_pack',&
              'release work storage')
  END IF

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  DEALLOCATE(label_atmname,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coord_pack','label_atmname')
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  !CALL stop_program("topology_coord_pack","checkpoint end")

END SUBROUTINE topology_coord_pack


!!*****
!******************************************************************************
!!****** topology/topology_pack [1.0] *
!!
!!   NAME
!!     topology_pack
!!
!!   SYNOPSIS
!!     Subroutine topology_pack(bond_a, bond_b, theta_a, theta_b, theta_c,&
!!         phi_a, phi_b, phi_c, phi_d, label_atmname, label_resname,&
!!         atm_mass, atm_charge, map_mol_typ, map_mol_num, mol_setup, setup,&
!!         globenv)
!!       Implicit None
!!       Type(molecule_type), Dimension(:), Pointer:: mol_setup
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
!!     End Subroutine topology_pack
!!
!!   FUNCTION
!!     topology pack
!!
!!   AUTHOR
!!     IKUO 08.17.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE topology_pack(bond_a,bond_b,theta_a,theta_b,theta_c,&
                     phi_a,phi_b,phi_c,phi_d,&
                     label_atmname,label_resname,&
                     atm_mass,atm_charge,&
                     map_mol_typ,map_mol_num,&
                     mol_setup,setup,globenv)

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
  CHARACTER (LEN=20), POINTER :: label_resname(:)
  CHARACTER (LEN=20), POINTER :: label_atmname(:)
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
  natom  = SIZE(map_mol_typ)
  nbond  = SIZE(bond_a)
  ntheta = SIZE(theta_a)
  nphi   = SIZE(phi_a)
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  WRITE(iw,*) 'Running the topology_pack'
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  nmol_typ = map_mol_typ(natom)
  IF(nmol_typ < 0) THEN
    CALL stop_program ( "topology_pack", "no moleculular type defined")
  ELSE
    !Allocate space for the molecules
    ALLOCATE(mol_setup(nmol_typ),STAT=ierror)
    IF(ierror/=0) CALL stop_memory("topology_pack","mol_setup",nmol_typ)
  END IF
!------------------------------------------------------------------------------
  DO i=1,nmol_typ
    CALL init_molecule_type(mol_setup(i))
    CALL find_boundary(map_mol_typ,natom,first,last,i)
    mol_setup(i)%molpar%name = label_resname(first)
    mol_setup(i)%num_mol=MAXVAL(map_mol_num(first:last))

    !--------------------------------------------------------------------------
    !set default nhcopt
    mol_setup(i)%nhcopt = "GLOBAL"
    !--------------------------------------------------------------------------
    !Need to be generalized!!!
    !Figure out the number of atoms and atom types-----------------------------
    j = last - first + 1
    j = j / map_mol_num(last)
    mol_setup(i)%molpar%natom=j
    mol_setup(i)%molpar%natom_type=1
    DO j=1,mol_setup(i)%molpar%natom-1
      IF(.NOT.(label_atmname(first+j-1)==label_atmname(first+j))) THEN
        mol_setup(i)%molpar%natom_type = mol_setup(i)%molpar%natom_type + 1
      END IF
    END DO
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !get the memory first
    ALLOCATE(mol_setup(i)%molpar%aname(mol_setup(i)%molpar%natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('topology_pack','aname',mol_setup(i)%molpar%natom)
    ALLOCATE(mol_setup(i)%molpar%alabel(mol_setup(i)%molpar%natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('topology_pack','alabel',mol_setup(i)%molpar%natom)
    ALLOCATE(mol_setup(i)%molpar%aweight(mol_setup(i)%molpar%natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('topology_pack','aweight',mol_setup(i)%molpar%natom)
    ALLOCATE(mol_setup(i)%molpar%acharge(mol_setup(i)%molpar%natom),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('topology_pack','acharge',mol_setup(i)%molpar%natom)
    DO j=1,mol_setup(i)%molpar%natom
      mol_setup(i)%molpar%aname(j)   = ''
      mol_setup(i)%molpar%aweight(j) = 0.0
      mol_setup(i)%molpar%acharge(j) = 0.0
      !Stuff the aname,alabel,aweight,acharge
      mol_setup(i)%molpar%aname(j)   = label_atmname(first+j-1)
      mol_setup(i)%molpar%aweight(j) = atm_mass(first+j-1)
      mol_setup(i)%molpar%acharge(j) = atm_charge(first+j-1)
    END DO
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Stuff the bond info
    ALLOCATE(checkme(nbond),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('topology_pack','checkme',nbond)
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
      ('topology_pack','mol_setup(i)%molpar%bonds',nbond)
    k=0
    DO j=1,nbond
      IF (checkme(j)==1) THEN
        k=k+1
        !-------------------------------------------------------------------
        !Check and make sure
        mol_setup(i)%molpar%bonds(1,k) = bond_a(j)
        mol_setup(i)%molpar%bonds(2,k) = bond_b(j)
        mol_setup(i)%molpar%bonds(3,k) = j
        !-------------------------------------------------------------------
      END IF
    END DO
    IF (k /= mol_setup(i)%molpar%nbonds) THEN
      !ERROR STOP PROGRAM
      CALL stop_program ( "topology_pack", "nbonds counting error")
    END IF
    DEALLOCATE(checkme,STAT=ierror)
    IF (ierror /= 0) CALL stop_memory ('read_molecule_section','checkme')
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Stuff the theata info
    ALLOCATE(checkme(ntheta),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('topology_pack','checkme',ntheta)
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
      ('topology_pack','mol_setup(i)%molpar%bends',&
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
        mol_setup(i)%molpar%bends(4,k) = j
        !-------------------------------------------------------------------
      END IF
    END DO
    IF (k /= mol_setup(i)%molpar%nbends) THEN
      !ERROR STOP PROGRAM
      CALL stop_program ( "topology_pack", "nbends counting error")
    END IF
    DEALLOCATE(checkme,STAT=ierror)
    IF (ierror /= 0) CALL stop_memory ('read_molecule_section','checkme')
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Stuff the torsion info
    ALLOCATE(checkme(nphi),STAT=ierror)
    IF(ierror/=0) CALL stop_memory &
      ('topology_pack','checkme',nphi)
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
        mol_setup(i)%molpar%torsions(5,k) = j
        !-------------------------------------------------------------------
      END IF
    END DO
    IF (k /= mol_setup(i)%molpar%ntorsions) THEN
      !ERROR STOP PROGRAM
      CALL stop_program ( "topology_pack", "ntorsions counting error")
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

  DEALLOCATE(label_resname,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','atm1')
  DEALLOCATE(atm_charge,atm_mass,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','atm2')
  DEALLOCATE(map_mol_typ,map_mol_num,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','map_mol_')
  DEALLOCATE(bond_a,bond_b,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','bond_')
  DEALLOCATE(theta_a,theta_b,theta_c,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','theta_')
  DEALLOCATE(phi_a,phi_b,phi_c,phi_d,STAT=ierror)
  IF ( ierror /= 0 ) CALL stop_memory ('read_molecule_section','phi_')

END SUBROUTINE topology_pack

END MODULE topology_util
