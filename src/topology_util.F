!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
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
!!    jgh (23-05-2004) Last atom of molecule information added
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
                                             allocate_potential,&
                                             set_potential
  USE force_field_types,               ONLY: allocate_bond_kind_set,&
                                             allocate_bend_kind_set,&
                                             allocate_torsion_kind_set,&
                                             bend_kind_type,&
                                             bond_kind_type,&
                                             torsion_kind_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl,&
                                             default_string_length
  USE memory_utilities,                ONLY: reallocate
  USE molecule_kind_types,             ONLY: allocate_molecule_kind_set,&
                                             atom_type,&
                                             bend_type,&
                                             bond_type,&
                                             dist_constraint_type,&
                                             g3x3_constraint_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type,&
                                             set_molecule_kind,&
                                             torsion_type
  USE molecule_types_new,              ONLY: allocate_molecule_set,&
                                             get_molecule,&
                                             local_constraint_type,&
                                             local_dist_constraint_type,&
                                             local_g3x3_constraint_type,&
                                             local_molecule_type,&
                                             molecule_type,&
                                             set_molecule,&
                                             set_molecule_set
  USE nrutil,                          ONLY: swap
  USE pair_potential_types,            ONLY: pair_potential_type
  USE particle_types,                  ONLY: allocate_particle_set,&
                                             particle_type
  USE periodic_table,                  ONLY: ptable
  USE physcon,                         ONLY: massunit
  USE simulation_cell,                 ONLY: pbc
  USE string_utilities,                ONLY: integer_to_string,&
                                             str_search
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE topology_types,                  ONLY: topology_parameters_type,&
                                             atom_info_type,&
                                             connectivity_info_type,&
                                             constraint_info_type,&
                                             force_field_amber_info_type,&
                                             force_field_charmm_info_type,&
                                             force_field_input_info_type

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "topology_util"


  PRIVATE
  PUBLIC :: find_boundary,&
            find_boundary2,&
            find_boundary3,&
            sort_bond,&
            sort_bend,&
            sort_dihe,&
            topology_reorder_atoms,&
            topology_connectivity_pack,&
            topology_constraint_pack,&
            topology_coordinate_pack,&
            topology_force_field_pack,&
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
!!     topology_reorder_atoms
!!
!!
!!   FUNCTION
!!
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE topology_reorder_atoms(topology,globenv)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_reorder_atoms', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: PM1, PM2, PM3, PM4

  TYPE(atom_info_type),POINTER     :: atom_info
  TYPE(connectivity_info_type),POINTER     :: conn_info

    CHARACTER (LEN=default_string_length), POINTER   :: tlabel_resname(:)
    CHARACTER (LEN=default_string_length), POINTER   :: tlabel_atmname(:)
    REAL (dbl), DIMENSION(:,:), POINTER :: tr
    INTEGER, POINTER :: tmap_mol_typ(:)
    INTEGER, POINTER :: tmap_mol_num(:)
    REAL (dbl), POINTER :: tatm_charge(:)
    REAL (dbl), POINTER :: tatm_mass(:)

    INTEGER      :: istat,iw
    INTEGER      :: location,iatm_count,imol_count
    INTEGER      :: ityp,inum,iatm,ntyp,nnum,natom
    INTEGER, POINTER :: new_position(:)
    INTEGER, POINTER :: mnum(:)
    INTEGER, POINTER :: mtyp(:)

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

  IF(PM1)  WRITE(iw,*) "    Entering topology_reorder_atoms"

  atom_info => topology%atom_info
  conn_info => topology%conn_info

  natom = topology%natoms

  NULLIFY(mtyp,mnum,new_position)
  NULLIFY(tlabel_resname,tlabel_atmname)
  NULLIFY(tr,tatm_charge,tatm_mass)
  NULLIFY(tmap_mol_typ,tmap_mol_num)

  !ALLOCATE all the temporary arrays needed for this routine
  ALLOCATE(mtyp(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','mtyp')
  ALLOCATE(mnum(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','mnum')
  ALLOCATE(new_position(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','new_position')
  ALLOCATE(tlabel_resname(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','tlabel_resname')
  ALLOCATE(tlabel_atmname(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','tlabel_atmname')
  ALLOCATE(tr(3,natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','tr')
  ALLOCATE(tatm_charge(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','tatm_charge')
  ALLOCATE(tatm_mass(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','tatm_mass')
  ALLOCATE(tmap_mol_typ(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','tmap_mol_typ')
  ALLOCATE(tmap_mol_num(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','tmap_mol_num')

  DO iatm=1,natom
    mtyp(iatm) = atom_info%map_mol_typ(iatm)
    mnum(iatm) = atom_info%map_mol_num(iatm)
  END DO
  ntyp = MAXVAL(mtyp)
  nnum = MAXVAL(mnum)

  imol_count = 0
  iatm_count = 0

  !Let's determine the new order that we want everything to be in
  DO ityp=1,ntyp
    DO inum=1,nnum
      imol_count = imol_count + 1
      DO iatm=1,natom
        IF ( (ityp==mtyp(iatm)) .AND. (inum==mnum(iatm)) ) THEN
          iatm_count = iatm_count + 1
          new_position(iatm_count) = iatm
          CALL integer_to_string(imol_count,atom_info%label_resname(iatm))
          atom_info%label_resname(iatm) = "R"//atom_info%label_resname(iatm)
        END IF
      END DO
    END DO
  END DO

  !Lets swap the atoms now
  DO iatm=1,natom
    location = new_position(iatm)
    tlabel_resname(iatm)       = atom_info%label_resname(location)
    tlabel_atmname(iatm)       = atom_info%label_atmname(location)
    tr(1,iatm)                 = atom_info%r(1,location)
    tr(2,iatm)                 = atom_info%r(2,location)
    tr(3,iatm)                 = atom_info%r(3,location)
    tatm_charge(iatm)          = atom_info%atm_charge(location)
    tatm_mass(iatm)            = atom_info%atm_mass(location)
    tmap_mol_typ(iatm)         = atom_info%map_mol_typ(location)
    tmap_mol_num(iatm)         = atom_info%map_mol_num(location)
  END DO
  DO iatm=1,natom
    atom_info%label_resname(iatm)   = tlabel_resname(iatm)
    atom_info%label_atmname(iatm)   = tlabel_atmname(iatm)
    atom_info%r(1,iatm)             = tr(1,iatm)
    atom_info%r(2,iatm)             = tr(2,iatm)
    atom_info%r(3,iatm)             = tr(3,iatm)
    atom_info%atm_charge(iatm)      = tatm_charge(iatm)
    atom_info%atm_mass(iatm)        = tatm_mass(iatm)
    atom_info%map_mol_typ(iatm)     = tmap_mol_typ(iatm)
    atom_info%map_mol_num(iatm)     = tmap_mol_num(iatm)
  END DO

  !DEALLOCATE all the temporary arrays needed for this routine
  DEALLOCATE(mtyp,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','mtyp')
  DEALLOCATE(mnum,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','mnum')
  DEALLOCATE(new_position,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','new_position')
  DEALLOCATE(tlabel_resname,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','tlabel_resname')
  DEALLOCATE(tlabel_atmname,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','tlabel_atmname')
  DEALLOCATE(tr,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','tr')
  DEALLOCATE(tatm_charge,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','tatm_charge')

  !DEALLOCATE all the stuff in topology structure so it can be re-runned
  DEALLOCATE(atom_info%map_mol_typ,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','map_mol_typ')
  DEALLOCATE(atom_info%map_mol_num,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','map_mol_num')
  DEALLOCATE(conn_info%bond_a,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','bond_a')
  DEALLOCATE(conn_info%bond_b,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','bond_b')
  DEALLOCATE(conn_info%theta_a,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','theta_a')
  DEALLOCATE(conn_info%theta_b,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','theta_b')
  DEALLOCATE(conn_info%theta_c,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','theta_c')
  DEALLOCATE(conn_info%phi_a,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','phi_a')
  DEALLOCATE(conn_info%phi_b,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','phi_b')
  DEALLOCATE(conn_info%phi_c,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','phi_c')
  DEALLOCATE(conn_info%phi_d,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_reorder_atoms','phi_d')


  IF(PM1) WRITE(iw,*) "    Exiting  topology_reorder_atoms"

  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("leaving "//routineN,globenv)
  END IF


!------------------------------------------------------------------------------
END SUBROUTINE topology_reorder_atoms

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

  TYPE(atom_info_type),POINTER     :: atom_info
  TYPE(connectivity_info_type),POINTER     :: conn_info

    INTEGER                                  :: first, iatom, iiatom, iii, &
                                                inum, istat, itype, iw, last, &
                                                natom, nbond, nnum, ntype
    INTEGER, POINTER                         :: checkme(:)
    LOGICAL                                  :: PM1, PM2, PM3, PM4

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

  atom_info => topology%atom_info
  conn_info => topology%conn_info

  ntype = 1
  itype = 0
  nnum = 1
  inum = 0

  natom = topology%natoms
  nbond = SIZE(conn_info%bond_a)
 
  ALLOCATE(checkme(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_generate_molecule','checkme')
  ALLOCATE (atom_info%map_mol_typ(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_generate_molecule','map_mol_typ')
  ALLOCATE (atom_info%map_mol_num(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_generate_molecule','map_mol_num')

  !Zero the arrays
  checkme(:) = 0
  atom_info%map_mol_typ(:) = 0
  atom_info%map_mol_num(:) = 0

  !Parse the atoms list to see how many different molecule types there are
  ntype = 1
  atom_info%map_mol_typ(1)=1
  DO iatom=2,natom
    IF(atom_info%label_molname(iatom-1)==atom_info%label_molname(iatom)) THEN
      atom_info%map_mol_typ(iatom)=ntype
    ELSE
      ntype=ntype+1
      atom_info%map_mol_typ(iatom)=ntype
    END IF
  END DO

  DO itype=1,ntype
    !search for first atom in the type list and last to cut down time
    CALL find_boundary(atom_info%map_mol_typ,natom,first,last,itype)
    IF(PM2) WRITE(iw,*) "      itype ",itype,first,last

    inum = 0

    DO iatom=first,last
      checkme(:) = 0

      IF(atom_info%map_mol_num(iatom)==0) THEN
        inum=inum+1
        checkme(iatom)=1
        DO WHILE (SUM(checkme)>0)
          DO iiatom=first,last
            IF((checkme(iiatom)==1)) THEN
              checkme(iiatom)=0
              atom_info%map_mol_num(iiatom)=inum
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
                IF((conn_info%bond_a(iii)==iiatom).AND.&
                  (atom_info%map_mol_num(conn_info%bond_b(iii))==0)) THEN 
                  !WRITE(*,*) "bond1 ",iiatom,conn_info%bond_b(iii)
                  checkme(conn_info%bond_b(iii))=1
                END IF
                IF((conn_info%bond_b(iii)==iiatom).AND.&
                  (atom_info%map_mol_num(conn_info%bond_a(iii))==0)) THEN
                  !WRITE(*,*) "bond2 ",iiatom,conn_info%bond_a(iii)
                  checkme(conn_info%bond_a(iii))=1
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
                 atom_info%map_mol_typ(iatom),atom_info%map_mol_num(iatom)
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
            CALL swap(phi_d(j),phi_d(k))
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
            CALL swap(phi_d(j),phi_d(k))
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
  TYPE(atom_info_type),POINTER     :: atom_info
  TYPE(connectivity_info_type),POINTER     :: conn_info

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

  atom_info => topology%atom_info
  conn_info => topology%conn_info

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
    atm_symbol(iatom) = ''
    IF(ASSOCIATED(atom_info%element)) THEN
      atm_symbol(iatom) = atom_info%element(iatom)
    END IF
    IF(atm_symbol(iatom) == '') THEN
      atm_symbol(iatom) = atom_info%label_atmname(iatom)
    END IF
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

        atom_info%atm_mass(iatom) = ptable(ielem)% amass
        !Need to figure out the correct charge!!!
        atom_info%atm_charge(iatom) = 0.0_dbl
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
      dr(1) = atom_info%r(1,iatom) - atom_info%r(1,jatom)
      dr(2) = atom_info%r(2,iatom) - atom_info%r(2,jatom)
      dr(3) = atom_info%r(3,iatom) - atom_info%r(3,jatom)

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
        IF(atom_info%label_molname(iatom)/=atom_info%label_molname(jatom)) THEN
          CALL stop_program ("topology_generate_bond",&
                             "bonds between different molecule types??? ")
        END IF

      END IF
    END DO
  END DO
!------------------------------------------------------------------------------
  !Copy the stuff from temporary array into the more permanent one
  ALLOCATE(conn_info%bond_a(nbond),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','bond_a',nbond)
  ALLOCATE(conn_info%bond_b(nbond),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bond','bond_b',nbond)
  DO ibond=1,nbond
    conn_info%bond_a(ibond) = tbond_a(ibond)
    conn_info%bond_b(ibond) = tbond_b(ibond)
  END DO
!------------------------------------------------------------------------------
  !Sort out the bond list
  DO ibond=1,nbond
    IF(conn_info%bond_a(ibond)>conn_info%bond_b(ibond)) THEN 
      CALL swap(conn_info%bond_a(ibond),conn_info%bond_b(ibond))
    END IF
  END DO
  !Need to update this to quicksort for systems greater then 20000 atoms
  DO i=1,nbond
    DO j=1,nbond
      IF(conn_info%bond_a(i)<conn_info%bond_a(j)) THEN
        CALL swap(conn_info%bond_a(i),conn_info%bond_a(j))
        CALL swap(conn_info%bond_b(i),conn_info%bond_b(j))
      END IF
    END DO
  END DO
  i=1
  DO WHILE (i<nbond)
    search = conn_info%bond_a(i)
    CALL find_boundary(conn_info%bond_a,nbond,first,last,search)
    IF(.NOT.(first==last)) THEN
      DO j=first,last
        DO k=first,last
          IF(conn_info%bond_b(j)<conn_info%bond_b(k)) THEN
            CALL swap(conn_info%bond_a(j),conn_info%bond_a(k))
            CALL swap(conn_info%bond_b(j),conn_info%bond_b(k))
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
    IF((conn_info%bond_a(ibond-1)==conn_info%bond_a(ibond)).AND.&
       (conn_info%bond_b(ibond-1)==conn_info%bond_b(ibond))) THEN
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
  TYPE(connectivity_info_type),POINTER     :: conn_info

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

  conn_info => topology%conn_info

  nbond = SIZE(conn_info%bond_a)
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
    search = conn_info%bond_a(ibond)
    CALL find_boundary(conn_info%bond_a,nbond,first,last,search)
    IF (first /= 0) THEN
      DO iibond=first,last
        IF (conn_info%bond_b(ibond)/=conn_info%bond_b(iibond)) THEN
          ntheta = ntheta + 1
          IF(ntheta > SIZE(tmp_a)) THEN
            CALL reallocate(tmp_a,1,ntheta+nblock)
            CALL reallocate(tmp_b,1,ntheta+nblock)
            CALL reallocate(tmp_c,1,ntheta+nblock)
          END IF
          tmp_a(ntheta) = conn_info%bond_b(ibond)
          tmp_b(ntheta) = search
          tmp_c(ntheta) = conn_info%bond_b(iibond)
        END IF
      END DO
    END IF
    !Case B
    search = conn_info%bond_b(ibond)
    CALL find_boundary(conn_info%bond_a,nbond,first,last,search)
    IF (first /= 0) THEN
      DO iibond=first,last
        IF(conn_info%bond_a(ibond)/=conn_info%bond_b(iibond)) THEN
          ntheta = ntheta + 1
          IF(ntheta > SIZE(tmp_a)) THEN
            CALL reallocate(tmp_a,1,ntheta+nblock)
            CALL reallocate(tmp_b,1,ntheta+nblock)
            CALL reallocate(tmp_c,1,ntheta+nblock)
          END IF
          tmp_a(ntheta) = conn_info%bond_a(ibond)
          tmp_b(ntheta) = search
          tmp_c(ntheta) = conn_info%bond_b(iibond)
        END IF
      END DO
    END IF
    !Case C
    search = conn_info%bond_b(ibond)
    CALL find_boundary(conn_info%bond_b,nbond,first,last,search)
    IF (first /= 0) THEN
      DO iibond=first,last
        IF (conn_info%bond_a(ibond)/=conn_info%bond_a(iibond)) THEN
          ntheta = ntheta + 1
          IF(ntheta > SIZE(tmp_a)) THEN
            CALL reallocate(tmp_a,1,ntheta+nblock)
            CALL reallocate(tmp_b,1,ntheta+nblock)
            CALL reallocate(tmp_c,1,ntheta+nblock)
          END IF
          tmp_a(ntheta) = conn_info%bond_a(ibond)
          tmp_b(ntheta) = search
          tmp_c(ntheta) = conn_info%bond_a(iibond)
        END IF
      END DO
    END IF
  END DO
!------------------------------------------------------------------------------
  !Sort all the bends first
  CALL sort_bend(ntheta,tmp_a,tmp_b,tmp_c)

  ALLOCATE(conn_info%theta_a(ntheta),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bend','theta_a',ntheta)
  ALLOCATE(conn_info%theta_b(ntheta),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bend','theta_b',ntheta)
  ALLOCATE(conn_info%theta_c(ntheta),STAT=istat)
  IF(istat/=0) CALL stop_memory('generate_bend','theta_c',ntheta)
  countme=1
  IF (ntheta > 0) THEN
    conn_info%theta_a(countme) = tmp_a(1)
    conn_info%theta_b(countme) = tmp_b(1)
    conn_info%theta_c(countme) = tmp_c(1)
    DO itheta=2,ntheta
      IF( .NOT.(tmp_a(itheta) == conn_info%theta_a(countme) .AND. &
          tmp_b(itheta) == conn_info%theta_b(countme) .AND. &
          tmp_c(itheta) == conn_info%theta_c(countme))) THEN
        countme = countme + 1
        conn_info%theta_a(countme) = tmp_a(itheta)
        conn_info%theta_b(countme) = tmp_b(itheta)
        conn_info%theta_c(countme) = tmp_c(itheta)
      END IF
    END DO
    ntheta = countme
    CALL REALLOCATE(conn_info%theta_a,1,ntheta)
    CALL REALLOCATE(conn_info%theta_b,1,ntheta)
    CALL REALLOCATE(conn_info%theta_c,1,ntheta)
  END IF
  

  !WRITE(*,*) "NTHETA = ",ntheta,SIZE(conn_info%theta_a)

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

  TYPE(connectivity_info_type),POINTER     :: conn_info

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

  conn_info => topology%conn_info


  ntheta = SIZE(conn_info%theta_a)
  nphi = 0
  !Stupid and inefficient but faster than alternative
  ALLOCATE(tmp_a(2*ntheta))
  ALLOCATE(tmp_b(2*ntheta))
  ALLOCATE(tmp_c(2*ntheta))
  ALLOCATE(tmp_d(2*ntheta))
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  ALLOCATE(conn_info%phi_a(nphi))
  ALLOCATE(conn_info%phi_b(nphi))
  ALLOCATE(conn_info%phi_c(nphi))
  ALLOCATE(conn_info%phi_d(nphi))
  DO iphi=1,nphi
    conn_info%phi_a(iphi) = tmp_a(iphi)
    conn_info%phi_b(iphi) = tmp_b(iphi)
    conn_info%phi_c(iphi) = tmp_c(iphi)
    conn_info%phi_d(iphi) = tmp_d(iphi)
  END DO
  nphi = SIZE(conn_info%phi_a)
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

  TYPE(atom_info_type),POINTER     :: atom_info
  TYPE(connectivity_info_type),POINTER     :: conn_info

    CHARACTER(LEN=default_string_length)     :: nhcopt
    INTEGER                                  :: atm_a, atm_b, atm_c, atm_d, &
                                                counter, first, i, istat, iw, &
                                                j, last, natom, nbond, nmol, &
                                                nphi, ntheta
    INTEGER, DIMENSION(:), POINTER           :: molecule_list
    INTEGER, DIMENSION(:), POINTER           :: checkme, first_list, last_list
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

  atom_info => topology%atom_info
  conn_info => topology%conn_info

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. Set the topology%[nmol_type,nmol]
  !-----------------------------------------------------------------------------
  natom  = topology%natoms
  topology%nmol_type = atom_info%map_mol_typ(natom)
  DO i=1,topology%nmol_type
    CALL find_boundary(atom_info%map_mol_typ,natom,first,last,i)
    topology%nmol=topology%nmol+atom_info%map_mol_num(last)
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
    CALL find_boundary(atom_info%map_mol_typ,natom,first,last,i)
    CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           kind_number=i,&
                           name=atom_info%label_molname(first),&
                           nhcopt=nhcopt)
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 5. Set the molecule_list for molecule_kind in molecule_kind_set
  !-----------------------------------------------------------------------------
  counter=0
  DO i=1,topology%nmol_type
    CALL find_boundary(atom_info%map_mol_typ,natom,first,last,i)
    nmol=atom_info%map_mol_num(last)
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
  ALLOCATE(last_list(topology%nmol),STAT=istat)
  IF(istat/=0) CALL stop_memory('coordinate_pack','last_list',topology%nmol)
  counter=0
  DO i=1,topology%nmol_type
    CALL find_boundary(atom_info%map_mol_typ,natom,first,last,i)
    nmol=atom_info%map_mol_num(last)
    DO j=1,nmol
      counter=counter+1
      CALL find_boundary2(atom_info%map_mol_typ,atom_info%map_mol_num,&
                          natom,first,last,i,j)
      first_list(counter)=first
      last_list(counter)=last
      IF(PM2) WRITE(iw,'(5X,A,4I5)') &
        "first_list=ikind,nmol,first_atom,last_atom",i,j,first,last
    END DO
  END DO
  CALL set_molecule_set(molecule_set,first_list,last_list)
  DEALLOCATE(first_list,STAT=istat)
  IF(istat/=0) CALL stop_memory('coordinate_pack','frist_list')
  DEALLOCATE(last_list,STAT=istat)
  IF(istat/=0) CALL stop_memory('coordinate_pack','last_list')


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
    CALL find_boundary(atom_info%map_mol_typ,topology%natoms,first,last,i)
    natom = (last-first+1)/ atom_info%map_mol_num(last)
    ALLOCATE(atom_list(natom),STAT=istat)
    IF(istat/=0) CALL stop_memory('connectivity_pack','atom_list',natom)
    DO j=1,natom
      !Atomic kind information will be filled in (PART 2)
      NULLIFY(atom_list(j)%atomic_kind)
      atom_list(j)%name=atom_info%label_atmname(j+first-1)
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
    CALL find_boundary(atom_info%map_mol_typ,topology%natoms,first,last,i)
    nmol=atom_info%map_mol_num(last)
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
    nbond  = SIZE(conn_info%bond_a)
    CALL find_boundary(atom_info%map_mol_typ,topology%natoms,first,last,i)
    ALLOCATE(checkme(nbond),STAT=istat)
    IF(istat/=0) CALL stop_memory('connectivity_pack','checkme',nbond)
    checkme(:)=0
    DO j=1,nbond
      atm_a = conn_info%bond_a(j)
      atm_b = conn_info%bond_b(j)
      IF( (atom_info%map_mol_typ(atm_a)==i).AND.&
          (atom_info%map_mol_num(atm_a)==1).AND.&
          (atom_info%map_mol_typ(atm_b)==i).AND.&
          (atom_info%map_mol_num(atm_b)==1) ) THEN
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
        bond_list(counter)%a=conn_info%bond_a(j)-first+1
        bond_list(counter)%b=conn_info%bond_b(j)-first+1
        !point this to the right bond_kind_type if using force field
        NULLIFY(bond_list(j)%bond_kind)
        IF(PM2) THEN
          WRITE(iw,'(7X,A,I3,1X,A,I5,I5,1X,A,I5,I5)') "molecule_kind",i,"bond",&
                      conn_info%bond_a(j),&
                      conn_info%bond_b(j),&
                      "offset number at",&
                      conn_info%bond_a(j)-first+1,&
                      conn_info%bond_b(j)-first+1
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
  IF(counter/=SIZE(conn_info%bond_a)) THEN
    IF(PM2) WRITE(iw,*) "    checksum vs nbond ",counter,SiZE(conn_info%bond_a)
    CALL stop_program ( "connectivity_pack", "nbond counting error")
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 11. Set the molecule_kind%[nbend,bend_list] via set_molecule_kind
  !-----------------------------------------------------------------------------
  DO i=1,topology%nmol_type
    molecule_kind => molecule_kind_set(i)
    ntheta  = SIZE(conn_info%theta_a)
    CALL find_boundary(atom_info%map_mol_typ,topology%natoms,first,last,i)
    ALLOCATE(checkme(ntheta),STAT=istat)
    IF(istat/=0) CALL stop_memory('connectivity_pack','checkme',ntheta)
    checkme(:)=0
    DO j=1,ntheta
      atm_a = conn_info%theta_a(j)
      atm_b = conn_info%theta_b(j)
      atm_c = conn_info%theta_c(j)
      IF( (atom_info%map_mol_typ(atm_a)==i).AND.&
          (atom_info%map_mol_num(atm_a)==1).AND.&
          (atom_info%map_mol_typ(atm_b)==i).AND.&
          (atom_info%map_mol_num(atm_b)==1).AND.&
          (atom_info%map_mol_typ(atm_c)==i).AND.&
          (atom_info%map_mol_num(atm_c)==1) ) THEN
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
        bend_list(counter)%a=conn_info%theta_a(j)-first+1
        bend_list(counter)%b=conn_info%theta_b(j)-first+1
        bend_list(counter)%c=conn_info%theta_c(j)-first+1
        !point this to the right bond_kind_type if using force field
        NULLIFY(bend_list(j)%bend_kind)
        IF(PM2) THEN
          WRITE(iw,'(7X,A,I3,1X,A,I5,I5,I5,1X,A,I5,I5,I5)') &
                      "molecule_kind",i,"bend",&
                      conn_info%theta_a(j),&
                      conn_info%theta_b(j),&
                      conn_info%theta_c(j),&
                      "offset number at",&
                      conn_info%theta_a(j)-first+1,&
                      conn_info%theta_b(j)-first+1,&
                      conn_info%theta_c(j)-first+1
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
  IF(counter/=SIZE(conn_info%theta_a)) THEN
    IF(PM2) WRITE(iw,*) "    checksum vs nbend ",counter,SiZE(conn_info%theta_a)
    CALL stop_program ( "connectivity_pack", "nbend counting error")
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 12. Set the molecule_kind%[ntorsion,torsion_list] via set_molecule_kind
  !-----------------------------------------------------------------------------
  DO i=1,topology%nmol_type
    molecule_kind => molecule_kind_set(i)
    nphi  = SIZE(conn_info%phi_a)
    CALL find_boundary(atom_info%map_mol_typ,topology%natoms,first,last,i)
    ALLOCATE(checkme(nphi),STAT=istat)
    IF(istat/=0) CALL stop_memory('connectivity_pack','checkme',nphi)
    checkme(:)=0
    DO j=1,nphi
      atm_a = conn_info%phi_a(j)
      atm_b = conn_info%phi_b(j)
      atm_c = conn_info%phi_c(j)
      atm_d = conn_info%phi_d(j)
      IF( (atom_info%map_mol_typ(atm_a)==i).AND.&
          (atom_info%map_mol_num(atm_a)==1).AND.&
          (atom_info%map_mol_typ(atm_b)==i).AND.&
          (atom_info%map_mol_num(atm_b)==1).AND.&
          (atom_info%map_mol_typ(atm_c)==i).AND.&
          (atom_info%map_mol_num(atm_c)==1).AND.&
          (atom_info%map_mol_typ(atm_d)==i).AND.&
          (atom_info%map_mol_num(atm_d)==1) ) THEN
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
        torsion_list(counter)%a=conn_info%phi_a(j)-first+1
        torsion_list(counter)%b=conn_info%phi_b(j)-first+1
        torsion_list(counter)%c=conn_info%phi_c(j)-first+1
        torsion_list(counter)%d=conn_info%phi_d(j)-first+1
        !point this to the right bond_kind_type if using force field
        NULLIFY(torsion_list(j)%torsion_kind)
        IF(PM2) THEN
          WRITE(iw,'(7X,A,I3,1X,A,I5,I5,I5,I5,1X,A,I5,I5,I5,I5)') &
                      "molecule_kind",i,"torsion",&
                      conn_info%phi_a(j),&
                      conn_info%phi_b(j),&
                      conn_info%phi_c(j),&
                      conn_info%phi_d(j),&
                      "offset number at",&
                      conn_info%phi_a(j)-first+1,&
                      conn_info%phi_b(j)-first+1,&
                      conn_info%phi_c(j)-first+1,&
                      conn_info%phi_d(j)-first+1
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
  IF(counter/=SIZE(conn_info%phi_a)) THEN
    IF(PM2) WRITE(iw,*) "    checksum vs nphi ",counter,SiZE(conn_info%phi_a)
    CALL stop_program ( "connectivity_pack", "nphi counting error")
  END IF



!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
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

  TYPE(atom_info_type),POINTER     :: atom_info

    CHARACTER(LEN=default_string_length), &
      POINTER                                :: work(:)
    INTEGER                                  :: counter, first, i, ikind, &
                                                istat, iw, j, k, last, natom
    INTEGER, DIMENSION(:), POINTER           :: atomkind,atomlist, kind_of, &
                                                natom_of_kind
    LOGICAL                                  :: PM1, PM2, PM3, PM4
    REAL(dbl), DIMENSION(:), POINTER         :: charge, mass
    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
    TYPE(bend_type), DIMENSION(:), POINTER   :: bend_list
    TYPE(bond_type), DIMENSION(:), POINTER   :: bond_list
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(elp_potential_type), POINTER        :: elp_potential
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(torsion_type), DIMENSION(:), &
      POINTER                                :: torsion_list
    INTEGER, DIMENSION(:), POINTER :: list
    INTEGER                                  :: imolecule,first_atom,last_atom,&
                                                nbond,nbend,ntorsion

    CHARACTER(LEN=default_string_length)     :: atmname
    LOGICAL                                  :: found


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


  atom_info => topology%atom_info

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. Determine topology%[natom_type,atom_names] and save mass(natom_type) 
  !-----------------------------------------------------------------------------
  counter = 0
  ALLOCATE(work(topology%natoms),mass(topology%natoms),STAT=istat)
  IF (istat/=0) CALL stop_memory ('coordinate_pack','work,mass',topology%natoms)
  work(:)=""
  DO i = 1, topology%nmol_type
    CALL find_boundary(atom_info%map_mol_typ,topology%natoms,first,last,i)
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           natom=natom,atom_list=atom_list)
    DO j = 1, natom
      IF(str_search(work,counter,atom_list(j)%name)==0)THEN
         counter = counter + 1
         work(counter) = atom_list(j)%name
         mass(counter) = atom_info%atm_mass(first+j-1)
        IF(PM2) WRITE(iw,'(7X,A,1X,A5,F10.5)') &
                "NEW ATOMIC KIND",work(counter),mass(counter)
      END IF
    END DO
  END DO
  topology%natom_type = counter
  ALLOCATE(atom_info%atom_names(topology%natom_type),STAT=istat)
  IF(istat/=0) CALL stop_memory ('coordinate_pack','atom_names')
  atom_info%atom_names(1:counter)=work(1:counter)
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
    mass(i) = mass(i) * massunit
    CALL set_atomic_kind(atomic_kind=atomic_kind,kind_number=i,&
                         name=atom_info%atom_names(i),mass=mass(i)) 
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
      IF(atom_info%atom_names(i)==atom_info%label_atmname(j)) THEN
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
  ! 7. Fill in some more infor for particle_set and lpi (local part info)
  !-----------------------------------------------------------------------------
  DO i=1,topology%natoms
    ikind = kind_of(i)
    particle_set(i)%atomic_kind => atomic_kind_set(ikind)
    particle_set(i)%iatom = i
    particle_set(i)%r(:) = atom_info%r(:,i)
    NULLIFY(particle_set(i)%lpi)
    ALLOCATE(particle_set(i)%lpi)
    NULLIFY(particle_set(i)%lpi%exclusion)
    ALLOCATE(particle_set(i)%lpi%exclusion(SIZE(atomic_kind_set)))
  END DO
  DEALLOCATE(kind_of,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','kind_of')


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 8. Fill in the particle_set%lpi%exclusion list
  !-----------------------------------------------------------------------------
  ALLOCATE(atomkind(SIZE(particle_set)),STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','atomkind')
  DO i=1,SIZE(particle_set)
    atomic_kind => particle_set(i)%atomic_kind
    CALL get_atomic_kind(atomic_kind,kind_number=j)
    atomkind(i) = j
  END DO
  DO imolecule=1,SIZE(molecule_set)
    molecule => molecule_set(imolecule)
    CALL get_molecule(molecule,first_atom=first_atom,last_atom=last_atom)

    natom = last_atom - first_atom + 1

    ALLOCATE(atomlist(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','atomlist')

    DO i=1,natom

      atomlist(:) = 1 
      atomlist(i) = 0    !exclude yourself???

      CALL get_molecule_kind(molecule_kind,nbond=nbond,&
                             bond_list=bond_list)
      DO j=1,nbond
        IF( (bond_list(j)%a==i).OR. &
            (bond_list(j)%b==i) ) THEN
          atomlist(bond_list(j)%a) = 0
          atomlist(bond_list(j)%b) = 0
        END IF
      END DO
      CALL get_molecule_kind(molecule_kind,nbend=nbend,&
                             bend_list=bend_list)
      DO j=1,nbend
        IF( (bend_list(j)%a==i).OR. &
            (bend_list(j)%b==i).OR. &
            (bend_list(j)%c==i) ) THEN
          atomlist(bend_list(j)%a) = 0
          atomlist(bend_list(j)%b) = 0
          atomlist(bend_list(j)%c) = 0
        END IF
      END DO
      CALL get_molecule_kind(molecule_kind,ntorsion=ntorsion,&
                             torsion_list=torsion_list)
      DO j=1,ntorsion
        IF( (torsion_list(j)%a==i).OR. &
            (torsion_list(j)%b==i).OR. &
            (torsion_list(j)%c==i).OR. &
            (torsion_list(j)%d==i) ) THEN
          atomlist(torsion_list(j)%a) = 0
          atomlist(torsion_list(j)%b) = 0
          atomlist(torsion_list(j)%c) = 0
          atomlist(torsion_list(j)%d) = 0
        END IF
      END DO

      DO j=1,SIZE(atomic_kind_set)
        counter = 0
        NULLIFY(list)
        DO k=1,natom
          IF((atomkind(k+first_atom-1)==j).AND.(atomlist(k)==0)) THEN
            counter = counter + 1
            CALL reallocate(list,1,counter)
            list(counter) = k + first_atom - 1
          END IF
        END DO
        particle_set(i+first_atom-1)%lpi%exclusion(j)%list => list
      END DO
    END DO
    DEALLOCATE(atomlist,STAT=istat)
  END DO
  DEALLOCATE(atomkind,STAT=istat)
  IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','atomkind')


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 9. Set the atom_list for molecule_kind in molecule_kind_set (PART 2)
  !-----------------------------------------------------------------------------
    DO i=1,topology%nmol_type
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           natom=natom,atom_list=atom_list)
    DO j=1,natom
      DO k=1,topology%natom_type
        IF(atom_list(j)%name == atom_info%atom_names(k)) THEN
          atom_list(j)%atomic_kind => atomic_kind_set(k)
        END IF
      END DO
    END DO
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 10. Set up atomic_kind_set()%elp_potentail%qeff
  !-----------------------------------------------------------------------------
  IF(topology%ff) THEN
    ALLOCATE(charge(SIZE(atomic_kind_set)),STAT=istat)
    IF (istat /= 0) CALL stop_memory('topology_coordinate_pack','charge')
    DO i=1,SIZE(atomic_kind_set)
      atomic_kind => atomic_kind_set(i)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                        name=atmname) 
  
      found = .FALSE.
  
      DO j=1,SIZE(topology%atom_info%label_atmname)
        IF(TRIM(topology%atom_info%label_atmname(j))==atmname) THEN
          charge(i) = topology%atom_info%atm_charge(j)
        END IF
      END DO
  
!     ! loop over params from amber
!     IF(ASSOCIATED(topology%ffam_info%ff_charge_atm)) THEN
!       DO j=1,SIZE(topology%ffam_info%ff_charge_atm)
!         IF(TRIM(topology%ffam_info%ff_charge_atm(j))==atmname) THEN
!           charge(i) = topology%ffam_info%ff_charge(j)
!           IF(found .AND. PM1) &
!             WRITE(iw,*) "    Multiple Charge declarations: overwriting"
!           found = .TRUE.
!         END IF
!       END DO
!     END IF
! 
!     ! loop over params from charmm
!     IF(ASSOCIATED(topology%ffch_info%ff_charge_atm)) THEN
!       DO j=1,SIZE(topology%ffch_info%ff_charge_atm)
!         IF(TRIM(topology%ffch_info%ff_charge_atm(j))==atmname) THEN
!           charge(i) = topology%ffch_info%ff_charge(j)
!           IF(found .AND. PM1) &
!             WRITE(iw,*) "    Multiple Charge declarations: overwriting"
!           found = .TRUE.
!         END IF
!       END DO
!     END IF
  
      ! always have the input param last to overwrite all the other ones
      IF(ASSOCIATED(topology%ffin_info%ff_charge_atm)) THEN
        DO j=1,SIZE(topology%ffin_info%ff_charge_atm)
          IF(TRIM(topology%ffin_info%ff_charge_atm(j))==atmname) THEN
            charge(i) = topology%ffin_info%ff_charge(j)
            IF(found .AND. PM1) &
              WRITE(iw,*) "    Multiple Charge declarations: overwriting"
            found = .TRUE.
          END IF
        END DO
      END IF

      IF(.NOT.found) THEN
        CALL stop_program("FF_coordinate_pack","missing qeff parm")
      END IF

      NULLIFY(elp_potential)
      CALL allocate_potential(elp_potential)
      CALL set_potential(potential=elp_potential,&
                         qeff=charge(i))
      CALL set_atomic_kind(atomic_kind=atomic_kind,&
                           elp_potential=elp_potential)
    END DO
    DEALLOCATE (charge,STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory('topology_coordinate_pack','release work storage')
    END IF
  END IF


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  IF(PM1) WRITE(iw,*) "  Exiting  topology_coordinate_pack"
  CALL write_checkpoint_information("leaving "//routineN,globenv)
  !CALL stop_program("topology_coordinate_pack","checkpoint end")

END SUBROUTINE topology_coordinate_pack


!******************************************************************************
!!****** topology_util/topology_constraint_pack [1.0] *
!!
!!   NAME
!!     topology_constraint_pack
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Pack in all the information needed for the constraints
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*** *************************************************************************

SUBROUTINE topology_constraint_pack ( molecule_kind_set,molecule_set,&
                                      topology,globenv)

  TYPE(molecule_kind_type), DIMENSION(:), &
    POINTER                                :: molecule_kind_set
  TYPE(molecule_type), DIMENSION(:), &
    POINTER                                :: molecule_set
  TYPE(topology_parameters_type), &
    INTENT(INOUT)                          :: topology
  TYPE(global_environment_type), &
    INTENT(IN), OPTIONAL                   :: globenv

  CHARACTER(len=*), PARAMETER :: routineN = 'topology_constraint_pack', &
    routineP = moduleN//':'//routineN

  TYPE(constraint_info_type),POINTER     :: cons_info

  INTEGER                                  :: istat, iw, i, j, k, natom
  INTEGER                                  :: counter,nmolecule,ndist,nbond
  INTEGER                                  :: ng3x3
  INTEGER                                  :: first, last,search1,search2
  LOGICAL                                  :: PM1, PM2, PM3, PM4

  INTEGER, DIMENSION(:), POINTER           :: molecule_list
  CHARACTER(LEN=default_string_length)  :: name,element_symbol
  TYPE(atomic_kind_type), POINTER        :: atomic_kind
  TYPE(molecule_kind_type), POINTER        :: molecule_kind
  TYPE(molecule_type), POINTER        :: molecule

  TYPE(dist_constraint_type), DIMENSION(:), POINTER :: dist_list
  TYPE(g3x3_constraint_type), DIMENSION(:), POINTER :: g3x3_list
  TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
  TYPE(bond_type), DIMENSION(:), POINTER            :: bond_list
  TYPE(local_constraint_type), POINTER       :: lci
  TYPE ( local_dist_constraint_type ), DIMENSION ( : ) , POINTER :: ldist
  TYPE ( local_g3x3_constraint_type ), DIMENSION ( : ) , POINTER :: lg3x3




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

  IF(PM1) WRITE(iw,*) "  Entering topology_constraint_pack"

  cons_info => topology%cons_info

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. NULLIFY the molecule_set(imol)%lci via set_molecule_set
  !-----------------------------------------------------------------------------
  DO i=1,topology%nmol
    molecule => molecule_set(i)
    NULLIFY(lci)
    ALLOCATE(lci,STAT=istat)
    IF(istat/=0) CALL stop_memory('constraint_pack','lci',1)
    NULLIFY(lci%ldist)
    NULLIFY(lci%lg3x3)
    CALL set_molecule(molecule,lci=lci)
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. Set info necessary for atom constraint
  !-----------------------------------------------------------------------------
  IF(topology%const_atom) THEN
    CALL stop_program("topology_constraint_pack","NO CONST_ATOM yet")
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 3. Add more stuff to cons_info%const_dist_* if constraint hydrogen is on
  !-----------------------------------------------------------------------------
  IF(topology%const_hydr) THEN
    ndist = cons_info%nconst_dist
    DO i=1,SIZE(molecule_kind_set)
      molecule_kind => molecule_kind_set(i)
      CALL get_molecule_kind(molecule_kind=molecule_kind,&
                             atom_list=atom_list,&
                             nbond=nbond,bond_list=bond_list)
      DO j=1,SIZE(atom_list)
        atomic_kind => atom_list(j)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name,element_symbol=element_symbol) 
        IF((TRIM(name)=="H").OR.(TRIM(element_symbol)=="H")) THEN
          DO k=1,nbond
            IF((bond_list(k)%a==j).OR.(bond_list(k)%b==j)) THEN
              ndist = ndist + 1
              CALL reallocate(cons_info%const_dist_mol,1,ndist)
              CALL reallocate(cons_info%const_dist_a,1,ndist)
              CALL reallocate(cons_info%const_dist_b,1,ndist)
              CALL reallocate(cons_info%const_dist_dab,1,ndist)
              cons_info%const_dist_mol(ndist)  = i
              cons_info%const_dist_a(ndist)    = bond_list(k)%a
              cons_info%const_dist_b(ndist)    = bond_list(k)%b
              cons_info%const_dist_dab(ndist)  = bond_list(k)%bond_kind%r0
            END IF
          END DO
        END IF
      END DO
    END DO
    cons_info%nconst_dist = ndist
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 4. sort and find duplicates cons_info%const_dist_*
  !-----------------------------------------------------------------------------
  IF(topology%const_dist) THEN
    ndist = cons_info%nconst_dist
    !Sort the first column
    DO i=1,ndist
      DO j=1,ndist
        IF(cons_info%const_dist_mol(i)<cons_info%const_dist_mol(j)) THEN
          CALL swap(cons_info%const_dist_mol(i),cons_info%const_dist_mol(j))
          CALL swap(cons_info%const_dist_a(i),cons_info%const_dist_a(j))
          CALL swap(cons_info%const_dist_b(i),cons_info%const_dist_b(j))
          CALL swap(cons_info%const_dist_dab(i),cons_info%const_dist_dab(j))
        END IF
      END DO
    END DO
    !Sort the second column
    i=1
    DO WHILE (i<ndist)
      search1 = cons_info%const_dist_mol(i)
      CALL find_boundary(cons_info%const_dist_mol,ndist,first,last,search1)
      IF(.NOT.(first==last)) THEN
        DO j=first,last
          DO k=first,last
            IF(cons_info%const_dist_a(j)<cons_info%const_dist_a(k)) THEN
              CALL swap(cons_info%const_dist_mol(j),cons_info%const_dist_mol(k))
              CALL swap(cons_info%const_dist_a(j),cons_info%const_dist_a(k))
              CALL swap(cons_info%const_dist_b(j),cons_info%const_dist_b(k))
              CALL swap(cons_info%const_dist_dab(j),cons_info%const_dist_dab(k))
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
    DO WHILE (i<ndist)
      search1 = cons_info%const_dist_mol(i)
      search2 = cons_info%const_dist_a(i)
      CALL find_boundary2(cons_info%const_dist_mol,cons_info%const_dist_a,&
                          ndist,first,last,search1,search2)
      IF(.NOT.(first==last)) THEN
        DO j=first,last
          DO k=first,last
            IF(cons_info%const_dist_b(j)<cons_info%const_dist_b(k)) THEN
              CALL swap(cons_info%const_dist_mol(j),cons_info%const_dist_mol(k))
              CALL swap(cons_info%const_dist_a(j),cons_info%const_dist_a(k))
              CALL swap(cons_info%const_dist_b(j),cons_info%const_dist_b(k))
              CALL swap(cons_info%const_dist_dab(j),cons_info%const_dist_dab(k))
            END IF
          END DO
        END DO
        i=last+1
      ELSE
        i=i+1
      END IF
    END DO
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 5. Set the distance constraint molecule_kind_set(ikind)%dist_list
  !-----------------------------------------------------------------------------
  IF(topology%const_dist) THEN
    DO i=1,SIZE(molecule_kind_set)
      molecule_kind => molecule_kind_set(i)
      CALL get_molecule_kind(molecule_kind=molecule_kind,&
                             nmolecule=nmolecule,molecule_list=molecule_list)
      CALL find_boundary(cons_info%const_dist_mol,cons_info%nconst_dist,&
                         first,last,i)
      ndist = last-first+1
      ALLOCATE(dist_list(ndist),STAT=istat)
      IF (istat /= 0) CALL stop_memory('topology_constraint_pack','dist_list')
      DO j=1,nmolecule
        molecule => molecule_set(molecule_list(j))
        ALLOCATE(ldist(ndist),STAT=istat)
        IF (istat /= 0) CALL stop_memory('topology_constraint_pack','ldist')
        DO k=1,ndist
          ldist(k)%del_lam = 0.0_dbl
          ldist(k)%lambda = 0.0_dbl
          ldist(k)%scale = 0.0_dbl
          ldist(k)%scale_old = 0.0_dbl
          ldist(k)%sigma = 0.0_dbl
          ldist(k)%f(:) = 0.0_dbl
          ldist(k)%ra_old(:) = 0.0_dbl
          ldist(k)%rb_old(:) = 0.0_dbl
        END DO
        CALL set_molecule(molecule=molecule,ldist=ldist)
      END DO
      ndist = 0
      DO j=first,last
        ndist = ndist + 1
        dist_list(ndist)%a   = cons_info%const_dist_a(j)
        dist_list(ndist)%b   = cons_info%const_dist_b(j)
        dist_list(ndist)%dab = cons_info%const_dist_dab(j)
      END DO
      CALL set_molecule_kind(molecule_kind,dist_list=dist_list,ndist=ndist)
    END DO
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 6. Set the group 3x3 constraint g3x3_list
  !-----------------------------------------------------------------------------
  IF(topology%const_33) THEN
    DO i=1,SIZE(molecule_kind_set)
      molecule_kind => molecule_kind_set(i)
      CALL get_molecule_kind(molecule_kind=molecule_kind,&
                             nmolecule=nmolecule,molecule_list=molecule_list)
      CALL find_boundary(cons_info%const_g33_mol,cons_info%nconst_g33,&
                         first,last,i)
      ng3x3 = last-first+1
      ALLOCATE(g3x3_list(ng3x3),STAT=istat)
      IF (istat /= 0) CALL stop_memory('topology_constraint_pack','g3x3_list')
      DO j=1,nmolecule
        molecule => molecule_set(molecule_list(j))
        ALLOCATE(lg3x3(ng3x3),STAT=istat)
        IF (istat /= 0) CALL stop_memory('topology_constraint_pack','lg3x3')
        DO k=1,ng3x3
          lg3x3(k)%scale = 0.0_dbl
          lg3x3(k)%scale_old = 0.0_dbl
          lg3x3(k)%fa = 0.0_dbl
          lg3x3(k)%fb = 0.0_dbl
          lg3x3(k)%fc = 0.0_dbl
          lg3x3(k)%ra_old = 0.0_dbl
          lg3x3(k)%rb_old = 0.0_dbl
          lg3x3(k)%rc_old = 0.0_dbl
          lg3x3(k)%va = 0.0_dbl
          lg3x3(k)%vb = 0.0_dbl
          lg3x3(k)%vc = 0.0_dbl
          lg3x3(k)%va_old = 0.0_dbl
          lg3x3(k)%vb_old = 0.0_dbl
          lg3x3(k)%vc_old = 0.0_dbl
          lg3x3(k)%del_lam = 0.0_dbl
          lg3x3(k)%lambda = 0.0_dbl
        END DO
        CALL set_molecule(molecule=molecule,lg3x3=lg3x3)
      END DO
      
      ng3x3 = 0
      DO j=1,SIZE(cons_info%const_g33_a)
        ng3x3 = ng3x3 + 1
        g3x3_list(ng3x3)%a   = cons_info%const_g33_a(j)
        g3x3_list(ng3x3)%b   = cons_info%const_g33_b(j)
        g3x3_list(ng3x3)%c   = cons_info%const_g33_c(j)
        g3x3_list(ng3x3)%dab = cons_info%const_g33_dab(j)
        g3x3_list(ng3x3)%dac = cons_info%const_g33_dac(j)
        g3x3_list(ng3x3)%dbc = cons_info%const_g33_dbc(j)
      END DO
      CALL set_molecule_kind(molecule_kind,ng3x3=ng3x3,g3x3_list=g3x3_list)
    END DO
  END IF


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  IF(PM1) WRITE(iw,*) "  Exiting  topology_constraint_pack"
  CALL write_checkpoint_information("leaving "//routineN,globenv)
  !CALL stop_program("topology_constraint_pack","checkpoint end")

END SUBROUTINE topology_constraint_pack




!******************************************************************************
!!****** topology_util/topology_force_field_pack [1.0] *
!!
!!   NAME
!!     topology_force_field_pack
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Pack in all the information needed for the force_field
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*** *************************************************************************

SUBROUTINE topology_force_field_pack (particle_set,atomic_kind_set,&
                                      molecule_kind_set,molecule_set,&
                                      topology,globenv,potparm)

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
  TYPE ( pair_potential_type ), &
    POINTER                                :: potparm ( :, : )

  CHARACTER(len=*), PARAMETER :: routineN = 'topology_force_field_pack', &
    routineP = moduleN//':'//routineN

  TYPE(force_field_input_info_type),POINTER     :: ffin_info
  TYPE(force_field_charmm_info_type),POINTER    :: ffch_info
  TYPE(force_field_amber_info_type),POINTER     :: ffam_info

  INTEGER                                  :: istat, iw, i, j, k
  LOGICAL                                  :: PM1, PM2, PM3, PM4

  CHARACTER(LEN=default_string_length)   :: name
  TYPE(atomic_kind_type), POINTER        :: atomic_kind
  TYPE(molecule_kind_type), POINTER      :: molecule_kind

  INTEGER                                  :: natom,nbond,nbend,ntorsion
  INTEGER                                  :: atm_a,atm_b,atm_c,atm_d
  CHARACTER(LEN=default_string_length)   :: name_atm_a,name_atm_b
  CHARACTER(LEN=default_string_length)   :: name_atm_c,name_atm_d
  CHARACTER(LEN=default_string_length)   :: name_atm_a2,name_atm_b2
  CHARACTER(LEN=default_string_length)   :: name_atm_c2,name_atm_d2
  TYPE(atom_type), DIMENSION(:), POINTER            :: atom_list
  TYPE(bond_type), DIMENSION(:), POINTER            :: bond_list
  TYPE(bend_type), DIMENSION(:), POINTER            :: bend_list
  TYPE(torsion_type), DIMENSION(:), POINTER         :: torsion_list
  LOGICAL                                  :: found
  INTEGER                                  :: counter
  INTEGER,POINTER                          :: map_bond_kind(:)
  INTEGER,POINTER                          :: map_bend_kind(:)
  INTEGER,POINTER                          :: map_torsion_kind(:)
  TYPE(bond_kind_type), DIMENSION(:), POINTER       :: bond_kind_set
  TYPE(bend_kind_type), DIMENSION(:), POINTER       :: bend_kind_set
  TYPE(torsion_kind_type), DIMENSION(:), POINTER    :: torsion_kind_set


  NULLIFY(atomic_kind,molecule_kind,atom_list,bond_list,bend_list,&
       torsion_list,map_bond_kind,map_bend_kind,map_torsion_kind,&
       bond_kind_set, bend_kind_set, torsion_kind_set)

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

  IF(PM1) WRITE(iw,*) "  Entering topology_force_field_pack"

  ffin_info => topology%ffin_info
  ffch_info => topology%ffch_info
  ffam_info => topology%ffam_info

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. Determine the number of unique bond kind and allocate bond_kind_set
  !-----------------------------------------------------------------------------
  DO i=1,SIZE(molecule_kind_set)
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           natom=natom,atom_list=atom_list,&
                           nbond=nbond,bond_list=bond_list)
    IF(nbond>0) THEN
      ALLOCATE(map_bond_kind(nbond),STAT=istat)
      IF (istat /= 0) CALL stop_memory('force_field_pack','map_bond_kind')
      counter=1
      map_bond_kind(1) = 1
      DO j=1,nbond
        atm_a = bond_list(j)%a
        atomic_kind => atom_list(atm_a)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_a)
        atm_b = bond_list(j)%b
        atomic_kind => atom_list(atm_b)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_b)
        found = .FALSE.
        DO k=1,j
          atm_a = bond_list(k)%a
          atomic_kind => atom_list(atm_a)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_a2)
          atm_b = bond_list(k)%b
          atomic_kind => atom_list(atm_b)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_b2)
          IF( ((TRIM(name_atm_a)==TRIM(name_atm_a2)) .AND. &
               (TRIM(name_atm_b)==TRIM(name_atm_b2))) .OR. &
              ((TRIM(name_atm_a)==TRIM(name_atm_b2)) .AND. &
               (TRIM(name_atm_b)==TRIM(name_atm_a2))) ) THEN
            found = .TRUE.
            map_bond_kind(j) = map_bond_kind(k)
            EXIT
          END IF
        END DO
        IF(.NOT.found) THEN
          counter=counter+1
          map_bond_kind(j) = counter
        END IF
      END DO
      CALL allocate_bond_kind_set(bond_kind_set,counter)
      DO j=1,nbond
        bond_list(j)%bond_kind => bond_kind_set(map_bond_kind(j))
      END DO
      CALL set_molecule_kind(molecule_kind=molecule_kind,&
                             bond_kind_set=bond_kind_set,bond_list=bond_list)
      DEALLOCATE(map_bond_kind,STAT=istat)
      IF (istat /= 0) CALL stop_memory('force_field_pack','map_bond_kind')
    END IF
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. Determine the number of unique bend kind and allocate bend_kind_set
  !-----------------------------------------------------------------------------
  DO i=1,SIZE(molecule_kind_set)
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           natom=natom,atom_list=atom_list,&
                           nbend=nbend,bend_list=bend_list)
    IF(nbend>0) THEN
      ALLOCATE(map_bend_kind(nbend),STAT=istat)
      IF (istat /= 0) CALL stop_memory('force_field_pack','map_bend_kind')
      counter=1
      map_bend_kind(1) = 1
      DO j=1,nbend
        atm_a = bend_list(j)%a
        atomic_kind => atom_list(atm_a)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_a)
        atm_b = bend_list(j)%b
        atomic_kind => atom_list(atm_b)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_b)
        atm_c = bend_list(j)%c
        atomic_kind => atom_list(atm_c)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_c)
        found = .FALSE.
        DO k=1,j
          atm_a = bend_list(k)%a
          atomic_kind => atom_list(atm_a)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_a2)
          atm_b = bend_list(k)%b
          atomic_kind => atom_list(atm_b)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_b2)
          atm_c = bend_list(k)%c
          atomic_kind => atom_list(atm_b)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_c2)
          IF( ((TRIM(name_atm_a)==TRIM(name_atm_a2)) .AND. &
               (TRIM(name_atm_b)==TRIM(name_atm_b2)) .AND. &
               (TRIM(name_atm_c)==TRIM(name_atm_c2))) .OR. &
              ((TRIM(name_atm_a)==TRIM(name_atm_c2)) .AND. &
               (TRIM(name_atm_b)==TRIM(name_atm_b2)) .AND. &
               (TRIM(name_atm_c)==TRIM(name_atm_a2))) ) THEN
            found = .TRUE.
            map_bend_kind(j) = map_bend_kind(k)
            EXIT
          END IF
        END DO
        IF(.NOT.found) THEN
          counter=counter+1
          map_bend_kind(j) = counter
        END IF
      END DO
      CALL allocate_bend_kind_set(bend_kind_set,counter)
      DO j=1,nbend
        bend_list(j)%bend_kind => bend_kind_set(map_bend_kind(j))
      END DO
      CALL set_molecule_kind(molecule_kind=molecule_kind,&
                             bend_kind_set=bend_kind_set,bend_list=bend_list)
      DEALLOCATE(map_bend_kind,STAT=istat)
      IF (istat /= 0) CALL stop_memory('force_field_pack','map_bend_kind')
    END IF
  END DO



  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 3. Determine the number of unique torsion kind and allocate torsion_kind_set
  !-----------------------------------------------------------------------------
  DO i=1,SIZE(molecule_kind_set)
    counter=0
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           natom=natom,atom_list=atom_list,&
                           ntorsion=ntorsion,torsion_list=torsion_list)
    IF(ntorsion>0) THEN
      ALLOCATE(map_torsion_kind(ntorsion),STAT=istat)
      IF (istat /= 0) CALL stop_memory('force_field_pack','map_torsion_kind')
      DO j=1,ntorsion
        atm_a = torsion_list(j)%a
        atomic_kind => atom_list(atm_a)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_a)
        atm_b = torsion_list(j)%b
        atomic_kind => atom_list(atm_b)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_b)
        atm_c = torsion_list(j)%c
        atomic_kind => atom_list(atm_c)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_c)
        atm_d = torsion_list(j)%d
        atomic_kind => atom_list(atm_d)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_d)
        found = .FALSE.
        DO k=1,j
          atm_a = torsion_list(k)%a
          atomic_kind => atom_list(atm_a)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_a2)
          atm_b = torsion_list(k)%b
          atomic_kind => atom_list(atm_b)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_b2)
          atm_c = torsion_list(k)%c
          atomic_kind => atom_list(atm_c)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_c2)
          atm_d = torsion_list(k)%d
          atomic_kind => atom_list(atm_d)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_d2)
          IF( ((TRIM(name_atm_a)==TRIM(name_atm_a2)) .AND. &
               (TRIM(name_atm_b)==TRIM(name_atm_b2)) .AND. &
               (TRIM(name_atm_c)==TRIM(name_atm_c2)) .AND. &
               (TRIM(name_atm_d)==TRIM(name_atm_d2))) .OR. &
              ((TRIM(name_atm_a)==TRIM(name_atm_d2)) .AND. &
               (TRIM(name_atm_b)==TRIM(name_atm_c2)) .AND. &
               (TRIM(name_atm_c)==TRIM(name_atm_b2)) .AND. &
               (TRIM(name_atm_d)==TRIM(name_atm_a2))) ) THEN
            found = .TRUE.
            map_torsion_kind(j) = map_torsion_kind(k)
            EXIT
          END IF
        END DO
        IF(.NOT.found) THEN
          counter=counter+1
          map_torsion_kind(j) = counter
        END IF
      END DO
      CALL allocate_torsion_kind_set(torsion_kind_set,counter)
      DO j=1,ntorsion
        torsion_list(j)%torsion_kind => torsion_kind_set(map_torsion_kind(j))
      END DO
      CALL set_molecule_kind(molecule_kind=molecule_kind,&
                             torsion_kind_set=torsion_kind_set,torsion_list=torsion_list)
      DEALLOCATE(map_torsion_kind,STAT=istat)
      IF (istat /= 0) CALL stop_memory('force_field_pack','map_torsion_kind')
    END IF
  END DO



  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 4. BONDS
  !-----------------------------------------------------------------------------
  DO i=1,SIZE(molecule_kind_set)
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           natom=natom,atom_list=atom_list,&
                           nbond=nbond,bond_list=bond_list)
    DO j=1,nbond
      atm_a = bond_list(j)%a
      atomic_kind => atom_list(atm_a)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_a)
      atm_b = bond_list(j)%b
      atomic_kind => atom_list(atm_b)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_b)
      found = .FALSE.

      ! loop over params from amber
!     IF(ASSOCIATED(ffam_info%ff_bond_a)) THEN
!       DO k=1,SIZE(ffam_info%ff_bond_a)
!         IF( ((TRIM(ffam_info%ff_bond_a(k))==TRIM(name_atm_a)) .AND. &
!              (TRIM(ffam_info%ff_bond_b(k))==TRIM(name_atm_b))) .OR. &
!             ((TRIM(ffam_info%ff_bond_a(k))==TRIM(name_atm_b)) .AND. &
!              (TRIM(ffam_info%ff_bond_b(k))==TRIM(name_atm_a))) ) THEN
!           bond_list(j)%bond_kind%k  = ffam_info%ff_bond_k(k)
!           bond_list(j)%bond_kind%r0 = ffam_info%ff_bond_r0(k)
!         IF(found .AND. PM1) &
!            WRITE(iw,*) "    Multiple Bond declarations: overwriting"
!           found = .TRUE.
!           EXIT
!         END IF
!       END DO
!     END IF

      ! loop over params from charmm
      IF(ASSOCIATED(ffch_info%ff_bond_a)) THEN
        DO k=1,SIZE(ffch_info%ff_bond_a)
          IF( ((TRIM(ffch_info%ff_bond_a(k))==TRIM(name_atm_a)) .AND. &
               (TRIM(ffch_info%ff_bond_b(k))==TRIM(name_atm_b))) .OR. &
              ((TRIM(ffch_info%ff_bond_a(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(ffch_info%ff_bond_b(k))==TRIM(name_atm_a))) ) THEN
            bond_list(j)%bond_kind%k  = ffch_info%ff_bond_k(k)
            bond_list(j)%bond_kind%r0 = ffch_info%ff_bond_r0(k)
          IF(found .AND. PM1) &
             WRITE(iw,*) "    Multiple Bond declarations: overwriting"
            found = .TRUE.
            EXIT
          END IF
        END DO
      END IF

      ! always have the input param last to overwrite all the other ones
      IF(ASSOCIATED(ffin_info%ff_bond_a)) THEN
        DO k=1,SIZE(ffin_info%ff_bond_a)
          IF( ((TRIM(ffin_info%ff_bond_a(k))==TRIM(name_atm_a)) .AND. &
               (TRIM(ffin_info%ff_bond_b(k))==TRIM(name_atm_b))) .OR. &
              ((TRIM(ffin_info%ff_bond_a(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(ffin_info%ff_bond_b(k))==TRIM(name_atm_a))) ) THEN
            bond_list(j)%bond_kind%k  = ffin_info%ff_bond_k(k)
            bond_list(j)%bond_kind%r0 = ffin_info%ff_bond_r0(k)
          IF(found .AND. PM1) &
             WRITE(iw,*) "    Multiple Bond declarations: overwriting"
            found = .TRUE.
            EXIT
          END IF
        END DO
      END IF

      IF(.NOT.found) CALL stop_program("FF_pack","missing bond parm")

    END DO
    CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           bond_list=bond_list)
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 5. BENDS
  !-----------------------------------------------------------------------------
  DO i=1,SIZE(molecule_kind_set)
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           natom=natom,atom_list=atom_list,&
                           nbend=nbend,bend_list=bend_list)
    DO j=1,nbend
      atm_a = bend_list(j)%a
      atomic_kind => atom_list(atm_a)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_a)
      atm_b = bend_list(j)%b
      atomic_kind => atom_list(atm_b)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_b)
      atm_c = bend_list(j)%c
      atomic_kind => atom_list(atm_c)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_c)
      found = .FALSE.

      ! loop over params from amber
!     IF(ASSOCIATED(ffam_info%ff_bend_a)) THEN
!       DO k=1,SIZE(ffam_info%ff_bend_a)
!         IF( ((TRIM(ffam_info%ff_bend_a(k))==TRIM(name_atm_a)) .AND. &
!              (TRIM(ffam_info%ff_bend_b(k))==TRIM(name_atm_b)) .AND. &
!              (TRIM(ffam_info%ff_bend_c(k))==TRIM(name_atm_c))) .OR. &
!             ((TRIM(ffam_info%ff_bend_a(k))==TRIM(name_atm_c)) .AND. &
!              (TRIM(ffam_info%ff_bend_b(k))==TRIM(name_atm_b)) .AND. &
!              (TRIM(ffam_info%ff_bend_c(k))==TRIM(name_atm_a))) ) THEN
!           bend_list(j)%bend_kind%k      = ffam_info%ff_bend_k(k)
!           bend_list(j)%bend_kind%theta0 = ffam_info%ff_bend_theta0(k)
!           IF(found .AND. PM1) &
!              WRITE(iw,*) "    Multiple Bend declarations: overwriting"
!           found = .TRUE.
!           EXIT
!         END IF
!       END DO
!     END IF

      ! loop over params from charmm
      IF(ASSOCIATED(ffch_info%ff_bend_a)) THEN
        DO k=1,SIZE(ffch_info%ff_bend_a)
          IF( ((TRIM(ffch_info%ff_bend_a(k))==TRIM(name_atm_a)) .AND. &
               (TRIM(ffch_info%ff_bend_b(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(ffch_info%ff_bend_c(k))==TRIM(name_atm_c))) .OR. &
              ((TRIM(ffch_info%ff_bend_a(k))==TRIM(name_atm_c)) .AND. &
               (TRIM(ffch_info%ff_bend_b(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(ffch_info%ff_bend_c(k))==TRIM(name_atm_a))) ) THEN
            bend_list(j)%bend_kind%k      = ffch_info%ff_bend_k(k)
            bend_list(j)%bend_kind%theta0 = ffch_info%ff_bend_theta0(k)
            IF(found .AND. PM1) &
               WRITE(iw,*) "    Multiple Bend declarations: overwriting"
            found = .TRUE.
            EXIT
          END IF
        END DO
      END IF

      ! always have the input param last to overwrite all the other ones
      IF(ASSOCIATED(ffin_info%ff_bend_a)) THEN
        DO k=1,SIZE(ffin_info%ff_bend_a)
          IF( ((TRIM(ffin_info%ff_bend_a(k))==TRIM(name_atm_a)) .AND. &
               (TRIM(ffin_info%ff_bend_b(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(ffin_info%ff_bend_c(k))==TRIM(name_atm_c))) .OR. &
              ((TRIM(ffin_info%ff_bend_a(k))==TRIM(name_atm_c)) .AND. &
               (TRIM(ffin_info%ff_bend_b(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(ffin_info%ff_bend_c(k))==TRIM(name_atm_a))) ) THEN
            bend_list(j)%bend_kind%k      = ffin_info%ff_bend_k(k)
            bend_list(j)%bend_kind%theta0 = ffin_info%ff_bend_theta0(k)
            IF(found .AND. PM1) &
               WRITE(iw,*) "    Multiple Bend declarations: overwriting"
            found = .TRUE.
            EXIT
          END IF
        END DO
      END IF

      IF(.NOT.found) CALL stop_program("FF_pack","missing bend parm")

    END DO
    CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           bend_list=bend_list)
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 6. TORSION
  !-----------------------------------------------------------------------------
  DO i=1,SIZE(molecule_kind_set)
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           natom=natom,atom_list=atom_list,&
                           ntorsion=ntorsion,torsion_list=torsion_list)
    DO j=1,ntorsion
      atm_a = torsion_list(j)%a
      atomic_kind => atom_list(atm_a)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_a)
      atm_b = torsion_list(j)%b
      atomic_kind => atom_list(atm_b)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_b)
      atm_c = torsion_list(j)%c
      atomic_kind => atom_list(atm_c)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_c)
      atm_d = torsion_list(j)%d
      atomic_kind => atom_list(atm_d)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_d)
      found = .FALSE.

      ! loop over params from amber
!     IF(ASSOCIATED(ffch_info%ff_torsion_a)) THEN
!       DO k=1,SIZE(ffch_info%ff_torsion_a)
!         IF( ((TRIM(ffch_info%ff_torsion_a(k))==TRIM(name_atm_a)) .AND. &
!              (TRIM(ffch_info%ff_torsion_b(k))==TRIM(name_atm_b)) .AND. &
!              (TRIM(ffch_info%ff_torsion_c(k))==TRIM(name_atm_c)) .AND. &
!              (TRIM(ffch_info%ff_torsion_d(k))==TRIM(name_atm_d))) .OR. &
!             ((TRIM(ffch_info%ff_torsion_a(k))==TRIM(name_atm_d)) .AND. &
!              (TRIM(ffch_info%ff_torsion_b(k))==TRIM(name_atm_c)) .AND. &
!              (TRIM(ffch_info%ff_torsion_c(k))==TRIM(name_atm_b)) .AND. &
!              (TRIM(ffch_info%ff_torsion_d(k))==TRIM(name_atm_a))) ) THEN
!           torsion_list(j)%torsion_kind%k      = ffch_info%ff_torsion_k(k)
!           torsion_list(j)%torsion_kind%phi0 = ffch_info%ff_torsion_phi0(k)
!           IF(found .AND. PM1) &
!              WRITE(iw,*) "    Multiple Torsion declarations: overwriting"
!           found = .TRUE.
!           EXIT
!         END IF
!       END DO
!     END IF

      ! loop over params from charmm
      IF(ASSOCIATED(ffch_info%ff_torsion_a)) THEN
        DO k=1,SIZE(ffch_info%ff_torsion_a)
          IF( ((TRIM(ffch_info%ff_torsion_a(k))==TRIM(name_atm_a)) .AND. &
               (TRIM(ffch_info%ff_torsion_b(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(ffch_info%ff_torsion_c(k))==TRIM(name_atm_c)) .AND. &
               (TRIM(ffch_info%ff_torsion_d(k))==TRIM(name_atm_d))) .OR. &
              ((TRIM(ffch_info%ff_torsion_a(k))==TRIM(name_atm_d)) .AND. &
               (TRIM(ffch_info%ff_torsion_b(k))==TRIM(name_atm_c)) .AND. &
               (TRIM(ffch_info%ff_torsion_c(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(ffch_info%ff_torsion_d(k))==TRIM(name_atm_a))) ) THEN
            torsion_list(j)%torsion_kind%k      = ffch_info%ff_torsion_k(k)
            torsion_list(j)%torsion_kind%phi0 = ffch_info%ff_torsion_phi0(k)
            IF(found .AND. PM1) &
               WRITE(iw,*) "    Multiple Torsion declarations: overwriting"
            found = .TRUE.
            EXIT
          END IF
        END DO
      END IF

      ! always have the input param last to overwrite all the other ones
      IF(ASSOCIATED(ffin_info%ff_torsion_a)) THEN
        DO k=1,SIZE(ffin_info%ff_torsion_a)
          IF( ((TRIM(ffin_info%ff_torsion_a(k))==TRIM(name_atm_a)) .AND. &
               (TRIM(ffin_info%ff_torsion_b(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(ffin_info%ff_torsion_c(k))==TRIM(name_atm_c)) .AND. &
               (TRIM(ffin_info%ff_torsion_d(k))==TRIM(name_atm_d))) .OR. &
              ((TRIM(ffin_info%ff_torsion_a(k))==TRIM(name_atm_d)) .AND. &
               (TRIM(ffin_info%ff_torsion_b(k))==TRIM(name_atm_c)) .AND. &
               (TRIM(ffin_info%ff_torsion_c(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(ffin_info%ff_torsion_d(k))==TRIM(name_atm_a))) ) THEN
            torsion_list(j)%torsion_kind%k      = ffin_info%ff_torsion_k(k)
            torsion_list(j)%torsion_kind%phi0 = ffin_info%ff_torsion_phi0(k)
            IF(found .AND. PM1) &
               WRITE(iw,*) "    Multiple Torsion declarations: overwriting"
            found = .TRUE.
            EXIT
          END IF
        END DO
      END IF

      IF(.NOT.found) CALL stop_program("FF_pack","missing torsion parm")

    END DO
    CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           torsion_list=torsion_list)
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 7. PARAMETERS ON POT
  !-----------------------------------------------------------------------------
  ALLOCATE(potparm(SIZE(atomic_kind_set),SIZE(atomic_kind_set)),STAT=istat)
  IF (istat /= 0) CALL stop_memory('ff_pack','potparm')
  DO i=1,SIZE(atomic_kind_set)
    atomic_kind => atom_list(i)%atomic_kind
    CALL get_atomic_kind(atomic_kind=atomic_kind,name=name_atm_a)
    DO j=1,SIZE(atomic_kind_set)
      atomic_kind => atom_list(j)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,name=name_atm_b)
      found = .FALSE.
      DO k=1,SIZE(ffin_info%ff_nonbonded)
        IF(PM2) WRITE(iw,*) "    TESTING ",TRIM(name_atm_a),TRIM(name_atm_b),&
            " with ",TRIM(ffin_info%ff_nonbonded(k)%at1),&
            TRIM(ffin_info%ff_nonbonded(k)%at2)
        IF(((TRIM(name_atm_a)==TRIM(ffin_info%ff_nonbonded(k)%at1))  .AND. &
            (TRIM(name_atm_b)==TRIM(ffin_info%ff_nonbonded(k)%at2))) .OR.  &
           ((TRIM(name_atm_b)==TRIM(ffin_info%ff_nonbonded(k)%at1))  .AND. &
            (TRIM(name_atm_a)==TRIM(ffin_info%ff_nonbonded(k)%at2))) ) THEN
          potparm(i,j) = ffin_info%ff_nonbonded(k)
          IF(PM2) WRITE(iw,*) "    FOUND ",TRIM(name_atm_a)," ",TRIM(name_atm_b)
          found = .TRUE.
          EXIT
        END IF
      END DO
      IF(.NOT.found) THEN
        WRITE(*,*) "NOT FOUND ",TRIM(name_atm_a),TRIM(name_atm_b)
        CALL stop_program("topology_force_field_pack","missing pair pot")
      END IF
    END DO
  END DO



!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  IF(PM1) WRITE(iw,*) "  Exiting  topology_force_field_pack"
  CALL write_checkpoint_information("leaving "//routineN,globenv)

END SUBROUTINE topology_force_field_pack


END MODULE topology_util
