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

  USE kinds,                            ONLY : dbl
  USE nrutil,                           ONLY : swap
  USE periodic_table,                   ONLY : ptable
  USE termination,                      ONLY : stop_program, &
                                               stop_memory

  PRIVATE
  PUBLIC :: topology_generate_molecule,azero,find_boundary

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     topology_generate_molecule
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
  CHARACTER (LEN=20) :: str1,str2

!------------------------------------------------------------------------------

  ntype = 1
  itype = 0
  nnum = 1
  inum = 0

  ALLOCATE(checkme(natom))
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
    !WRITE(*,*) "iatom(",iatom,") = ",map_mol_typ(iatom)
  END DO

  WRITE(*,*) "NUMBER OF MOLECULAR TYPE TOTAL = ",ntype

  DO itype=1,ntype
    !search for first atom in the type list and last to cut down time
    CALL find_boundary(map_mol_typ,natom,first,last,itype)
    !WRITE(*,*) "itype",itype,"from",first,"to",last

    inum = 0

    DO iatom=first,last
      CALL azero(checkme,natom)

      IF(map_mol_num(iatom)==0) THEN
        !WRITE(*,*) "iatom =",iatom

        inum=inum+1
        !WRITE(*,*) "inum now = ",inum

        checkme(iatom)=1
       
        DO WHILE (SUM(checkme))
          DO iiatom=first,last
            IF((checkme(iiatom)==1)) THEN
              !WRITE(*,*) "iiatom = ",iiatom
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
          !WRITE(*,*) "SUM = ",SUM(checkme)
        END DO
      END IF

    END DO
    WRITE(*,*) "MOL_TYPE ",atm_mol_name(first),"CONTAINS",inum,"MOLECULES"
  END DO

! DO iatom=1,natom
!   WRITE(*,*) "iatom",iatom,"",map_mol_typ(iatom),map_mol_num(iatom)
! END DO

  DEALLOCATE(checkme)

END SUBROUTINE topology_generate_molecule

!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     find_boundary
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

END SUBROUTINE 

!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     azero
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
!!     topology_generate_bonds
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
SUBROUTINE topology_generate_bonds(natom,x,y,z,nbond,bond_a,bond_b)
  IMPLICIT NONE
!Arguments
  INTEGER :: natom
  INTEGER, POINTER :: x(:),y(:),z(:)
  INTEGER :: nbond
  INTEGER, POINTER :: bond_a(:),bond_b(:)

!Locals
  INTEGER :: i,j,k
  INTEGER :: iatom,jatom,ibond
  INTEGER :: first,last,search
  INTEGER :: ielem,nelem
  INTEGER, POINTER :: tbond_a(:),tbond_b(:),tbond_c(:)
  REAL :: dx,dy,dz,r2,rbond,rbond2
  INTEGER, POINTER :: pos_x(:),pos_y(:),pos_z(:)
  REAL, POINTER :: covalent_radius(:)
  CHARACTER (LEN = 2),POINTER :: atm_symbol(:)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  CALL stop_program ("topology_generate_bonds",&
                     "Needs to be debug")

  !Get the temporary array that we need
  ALLOCATE(tbond_a(2*natom))
  ALLOCATE(tbond_b(2*natom))
  ALLOCATE(tbond_c(2*natom))
  ALLOCATE(pos_x(natom))
  ALLOCATE(pos_y(natom))
  ALLOCATE(pos_z(natom))
  ALLOCATE(atm_symbol(natom))
  ALLOCATE(covalent_radius(natom))
!------------------------------------------------------------------------------
  !Set up the arrays that we need for this routine
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  CALL azero(tbond_a,2*natom)
  CALL azero(tbond_b,2*natom)
  CALL azero(tbond_c,2*natom)

  DO iatom=1,natom
    pos_x(iatom) = x(iatom)
    pos_y(iatom) = y(iatom)
    pos_z(iatom) = z(iatom)
    
    nelem = SIZE(ptable)
    DO ielem=1,nelem 
      IF(ptable(ielem)%symbol==atm_symbol(iatom)) THEN
        covalent_radius(iatom) = ptable(ielem)%covalent_radius
        EXIT   !!!! Correct??? XXXXXXXXXXXXXXXXXXXX
      END IF
    END DO
  END DO
!------------------------------------------------------------------------------
  !Find bonds according to periodic table
  nbond = 0
  DO iatom=1,natom
    DO jatom=iatom+1,natom
      !Calculate distance between the two atoms
      !Put in PBC XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !Put in PBC XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !Put in PBC XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !MUST PASS IN CELL INFO
      dx = pos_x(iatom) - pos_x(jatom)
      dy = pos_x(iatom) - pos_x(jatom)
      dz = pos_x(iatom) - pos_x(jatom)
      r2 = dx*dx + dy*dy + dz*dz
      rbond = covalent_radius(iatom)+covalent_radius(jatom)
      rbond2 = rbond*rbond
      !Test the distance to the sum of the covalent radius
      IF(r2 < rbond2) THEN
        nbond=nbond+1
        tbond_a(nbond)=iatom
        tbond_b(nbond)=jatom
      END IF
    END DO
  END DO
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
!------------------------------------------------------------------------------
  !Sort out the duplicates... (Shouldn't have any but just to check)
  DO ibond=2,nbond
    IF((tbond_a(ibond-1)==tbond_a(ibond)).AND.&
       (tbond_b(ibond-1)==tbond_b(ibond))) THEN
      CALL stop_program ("topology_generate_bonds",&
                         "duplicate bonds")
    END IF
  END DO
!------------------------------------------------------------------------------
  DEALLOCATE(tbond_a)
  DEALLOCATE(tbond_b)
  DEALLOCATE(tbond_c)
  DEALLOCATE(pos_x)
  DEALLOCATE(pos_y)
  DEALLOCATE(pos_z)
  DEALLOCATE(atm_symbol)
  DEALLOCATE(covalent_radius)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE topology_generate_bonds


!
!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     topology_generate_bends
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
SUBROUTINE topology_generate_bends(nbond,bond_a,bond_b,&
                                   ntheta,theta_a,theta_b,theta_c)
  IMPLICIT NONE

!Arguments
  INTEGER :: nbond
  INTEGER, POINTER :: bond_a(:),bond_b(:)
  INTEGER :: ntheta
  INTEGER, POINTER :: theta_a(:),theta_b(:),theta_c(:)

!Locals
  INTEGER :: first,last,search
  INTEGER :: ibond,iibond,itheta
  INTEGER, POINTER :: tmp_a(:),tmp_b(:),tmp_c(:)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  ntheta = 0
  !Stupid and inefficient but faster than alternative
  ALLOCATE(tmp_a(2*nbond))
  ALLOCATE(tmp_b(2*nbond))
  ALLOCATE(tmp_c(2*nbond))
  CALL stop_program("topology_generate_bends","Need to be tested")
!------------------------------------------------------------------------------
  DO ibond=1,nbond
    search = bond_b(ibond)
    CALL find_boundary(bond_a,nbond,first,last,search)
    DO iibond=first,last
      ntheta = ntheta + 1
      tmp_a(ntheta) = bond_a(ibond)
      tmp_b(ntheta) = bond_a(iibond)
      tmp_c(ntheta) = bond_b(iibond)
    END DO
  END DO
!------------------------------------------------------------------------------
  DO itheta=1,ntheta
    theta_a(itheta) = tmp_a(itheta)
    theta_b(itheta) = tmp_b(itheta)
    theta_c(itheta) = tmp_c(itheta)
  END DO
!------------------------------------------------------------------------------
  DEALLOCATE(tmp_a)
  DEALLOCATE(tmp_b)
  DEALLOCATE(tmp_c)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE topology_generate_bends
 
 
 
!
!******************************************************************************
!!****** topology_util [1.0] *
!!
!!   NAME
!!     topology_generate_dihe
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
SUBROUTINE topology_generate_dihe(nbond,bond_a,bond_b,&
                                   nphi,phi_a,phi_b,phi_c,phi_d)
  IMPLICIT NONE

!Arguments
  INTEGER :: nbond
  INTEGER, POINTER :: bond_a(:),bond_b(:)
  INTEGER :: nphi
  INTEGER, POINTER :: phi_a(:),phi_b(:),phi_c(:),phi_d(:)

!Locals
  INTEGER :: first,last,search
  INTEGER :: firstb,lastb,searchb
  INTEGER :: ibond,iibond,iiibond,iphi
  INTEGER, POINTER :: tmp_a(:),tmp_b(:),tmp_c(:),tmp_d(:)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  nphi = 0
  !Stupid and inefficient but faster than alternative
  ALLOCATE(tmp_a(2*nbond))
  ALLOCATE(tmp_b(2*nbond))
  ALLOCATE(tmp_c(2*nbond))
  ALLOCATE(tmp_d(2*nbond))
  CALL stop_program("topology_generate_dihe","Need to be tested")
!------------------------------------------------------------------------------
  DO ibond=1,nbond
    search = bond_b(ibond)
    CALL find_boundary(bond_a,nbond,first,last,search)
    DO iibond=first,last
      searchb = bond_b(iibond)
      CALL find_boundary(bond_a,nbond,firstb,lastb,searchb)
      DO iiibond=firstb,lastb
        nphi = nphi + 1
        tmp_a(nphi) = bond_a(ibond)
        tmp_b(nphi) = bond_a(iibond)
        tmp_c(nphi) = bond_a(iiibond)
        tmp_d(nphi) = bond_b(iiibond)
      END DO
    END DO
  END DO
!------------------------------------------------------------------------------
  DO iphi=1,nphi
    phi_a(iphi) = tmp_a(iphi)
    phi_b(iphi) = tmp_b(iphi)
    phi_c(iphi) = tmp_c(iphi)
    phi_d(iphi) = tmp_d(iphi)
  END DO
!------------------------------------------------------------------------------
  DEALLOCATE(tmp_a)
  DEALLOCATE(tmp_b)
  DEALLOCATE(tmp_c)
  DEALLOCATE(tmp_d)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE topology_generate_dihe
 
END MODULE topology_util

!******************************************************************************
