!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/linklists [1.0] *
!!
!!   NAME
!!     linklists
!!
!!   FUNCTION
!!     Generate linked lists for intramolecular dependencies
!!
!!   AUTHOR
!!     CJM & JGH
!!
!!   MODIFICATION HISTORY
!!     Variable names changed (DG) 05-Dec-2000
!!
!!   SOURCE
!******************************************************************************

MODULE linklists

  USE kinds, ONLY : dbl
  USE molecule_types, ONLY : molecule_structure_type, linklist_atoms, &
       bond_parameters_type, linklist_bonds, bend_parameters_type, &
       linklist_bends, torsion_parameters_type, linklist_torsions, &
       linklist_dist_constraints, linklist_g3x3_constraints,  &
       linklist_g4x6_constraints, linklist_exclusion, particle_node_type, &
       molecule_type
  USE termination, ONLY : stop_memory, stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: atom_link_list, exclusion, bond_link_list, bend_link_list
  PUBLIC :: torsion_link_list, distconst_link_list, g3x3_link_list, g4x6_link_list
  PUBLIC :: bonds, bends, torsions, dist_constraints, g3x3_constraints, g4x6_constraints

  LOGICAL :: bonds, bends, torsions, dist_constraints, g3x3_constraints, g4x6_constraints

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** linklists/atom_link_list [1.0] *
!!
!!   NAME
!!     atom_link_list
!!
!!   FUNCTION
!!     Generates linked lists for intramolecular dependencies
!!
!!   AUTHOR
!!     CJM & JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE atom_link_list ( molecule, pnode )

  IMPLICIT NONE

! Arguments
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( INOUT ) :: &
       molecule
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( INOUT ), TARGET :: &
       pnode

! Locals
  INTEGER :: imol, iatom, jatom
  TYPE ( linklist_atoms ), POINTER :: current_ll_atom
  INTEGER :: ios

!------------------------------------------------------------------------------

  ALLOCATE ( current_ll_atom, STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory( 'atom_link_list', 'current_ll_atom', 0 )
  jatom = 1
  DO imol = 1, size(molecule)
     molecule(imol) %ll_atoms => current_ll_atom
     DO iatom = 1, molecule(imol) %natoms_mol
        current_ll_atom%part => pnode(jatom)
        ALLOCATE ( current_ll_atom%next, STAT = ios )
        IF ( ios /= 0 ) &
             CALL stop_memory( 'atom_link_list', 'current_ll_atom%next', 0 )
        current_ll_atom => current_ll_atom%next
        jatom = jatom + 1
     END DO
  END DO

END SUBROUTINE atom_link_list

!!*****
!******************************************************************************
!!****** linklists/exclusion [1.0] *
!!
!!   NAME
!!     exclusion
!!
!!   FUNCTION
!!     Generates linked lists for intramolecular dependencies
!!
!!   AUTHOR
!!     CJM & JGH
!!
!!   MODIFICATION HISTORY
!!     Debugged parallel code (JGH) 25-Nov-2000
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE exclusion(molecule,pnode)

! Routine to compute the exclusion list

  IMPLICIT NONE

! Arguments
  TYPE (molecule_structure_type ), INTENT ( IN ), DIMENSION ( : ) :: molecule
  TYPE (particle_node_type ), INTENT ( INOUT ), DIMENSION ( : ) :: pnode

! Locals
  TYPE (linklist_exclusion), POINTER :: current_excl
  TYPE (linklist_atoms), POINTER :: ll_atom
  TYPE (linklist_bonds), POINTER :: ll_bond
  TYPE (linklist_bends), POINTER :: ll_bend
  TYPE (linklist_torsions), POINTER :: ll_torsion
  TYPE (linklist_dist_constraints), POINTER :: ll_dist
  TYPE (linklist_g3x3_constraints), POINTER :: ll_g3x3
  INTEGER :: imol, iatom, jatom, ibond, ibend, itorsion, idist, i, ii
  INTEGER :: jnode, ig3x3, ios
  INTEGER, DIMENSION ( : ), ALLOCATABLE :: bond, bend, torsion
  INTEGER :: natoms_mol, nbonds_mol, nbends_mol, ntorsions_mol, &
             ndist_mol, n3x3_mol
  LOGICAL :: match

!------------------------------------------------------------------------------

! allocating local variables containting intra atoms
  ALLOCATE (bond(2),STAT=ios)
  IF ( ios /= 0 ) CALL stop_memory( 'exclusion','bond', 2 )
  ALLOCATE (bend(3),STAT=ios)
  IF ( ios /= 0 ) CALL stop_memory( 'exclusion','bend', 3 )
  ALLOCATE (torsion(4),STAT=ios)
  IF ( ios /= 0 ) CALL stop_memory( 'exclusion','torsion', 4 )

! allocating the exclusion pointer
  ALLOCATE (current_excl,STAT=ios)
  IF ( ios /= 0 ) CALL stop_memory( 'exclusion', 'current_excl', 0 )

! initialize the nexcl counter for all atoms
  pnode%nexcl = 0

! looping through all of the molecules
  jnode = 0

  MLOOP: DO imol = 1, size(molecule)

! assigning local variables
     natoms_mol = molecule(imol) %natoms_mol
     nbonds_mol = molecule(imol) %nbonds_mol
     nbends_mol = molecule(imol) %nbends_mol
     ntorsions_mol = molecule(imol) %ntorsions_mol
     ndist_mol = molecule(imol) %ndistconst_mol
     n3x3_mol = molecule(imol) %ng3x3const_mol

! pointing local pointers to the start of the list
     ll_atom => molecule(imol) %ll_atoms

! looping over all the atoms in molecule imol
     ALOOP: DO iatom = 1, natoms_mol
        jnode = jnode + 1

! finding the current atom in the molecule
        jatom = ll_atom%part%p%iatom
        pnode(jnode) %ex => current_excl
        ll_bond => molecule(imol) %ll_bonds

! doing the bonds
        BONDL: DO ibond = 1, nbonds_mol
           bond(1) = ll_bond%p1%iatom
           bond(2) = ll_bond%p2%iatom

! checking to see if jatom is bonded
           match=.FALSE.
           DO ii = 1, 2
              IF (jatom==bond(ii)) THEN
                 match = .TRUE.
                 EXIT
              END IF
           END DO
           BONDIF: IF (match) THEN
              DO ii = 1, 2
                 CALL checklist(pnode(jnode),match,jatom,bond(ii))
                 IF ( .NOT. match) THEN
                    i = getnode(molecule,bond(ii))
                    current_excl%p => pnode(i) %p
                    pnode(jnode) %nexcl = pnode(jnode) %nexcl + 1
                    ALLOCATE (current_excl%next,STAT=ios)
                    IF ( ios /= 0 ) CALL stop_memory ( 'exclusion', &
                         'current_excl%next', 0 )
                    current_excl => current_excl%next
                 END IF
              END DO
           END IF BONDIF
           ll_bond => ll_bond%next
        END DO BONDL

! doing the bends
        ll_bend => molecule(imol) %ll_bends
        BENDL: DO ibend = 1, nbends_mol
           bend(1) = ll_bend%p1%iatom
           bend(2) = ll_bend%p2%iatom
           bend(3) = ll_bend%p3%iatom
! checking to see if jatom is part of a bend
           DO ii = 1, 3
              IF (jatom==bend(ii)) THEN
                 match = .TRUE.
                 EXIT
              END IF
           END DO
           BENDIF: IF (match) THEN
              DO ii = 1, 3
                 CALL checklist(pnode(jnode),match,jatom,bend(ii))
                 IF ( .NOT. match) THEN
                    i = getnode(molecule,bend(ii))
                    current_excl%p => pnode(i) %p
                    pnode(jnode) %nexcl = pnode(jnode) %nexcl + 1
                    ALLOCATE (current_excl%next,STAT=ios)
                    IF ( ios /= 0 ) CALL stop_memory ( 'exclusion', &
                         'current_excl%next', 0 )
                    current_excl => current_excl%next
                 END IF
              END DO
           END IF BENDIF
           ll_bend => ll_bend%next
        END DO BENDL

! doing the torsions
        ll_torsion => molecule(imol) %ll_torsions
        TORSL: DO itorsion = 1, ntorsions_mol
           torsion(1) = ll_torsion%p1%iatom
           torsion(2) = ll_torsion%p2%iatom
           torsion(3) = ll_torsion%p3%iatom
           torsion(4) = ll_torsion%p4%iatom
! checking to see if jatom is part of a torsion
           DO ii = 1, 4
              IF (jatom==torsion(ii)) THEN
                 match = .TRUE.
                 EXIT
              END IF
           END DO
           TORSIF: IF (match) THEN
              DO ii = 1, 4
                 CALL checklist(pnode(jnode),match,jatom,torsion(ii))
                 IF ( .NOT. match) THEN
                    i = getnode(molecule,torsion(ii))
                    current_excl%p => pnode(i) %p
                    pnode(jnode) %nexcl = pnode(jnode) %nexcl + 1
                    ALLOCATE (current_excl%next,STAT=ios)
                    IF ( ios /= 0 ) CALL stop_memory ( 'exclusion', &
                         'current_excl%next', 0 )
                    current_excl => current_excl%next
                 END IF
              END DO
           END IF TORSIF
           ll_torsion => ll_torsion%next
        END DO TORSL

! doing the distance constraints
        ll_dist => molecule(imol) %ll_dist_const
        DISTL: DO idist = 1, ndist_mol
           bond(1) = ll_dist%p1%iatom
           bond(2) = ll_dist%p2%iatom

! checking to see if jatom is bonded
           DO ii = 1, 2
              IF (jatom==bond(ii)) THEN
                 match = .TRUE.
                 EXIT
              END IF
           END DO
           DISTIF: IF (match) THEN
              DO ii = 1, 2
                 CALL checklist(pnode(jnode),match,jatom,bond(ii))
                 IF ( .NOT. match) THEN
                    i = getnode(molecule,bond(ii))
                    current_excl%p => pnode(i) %p
                    pnode(jnode) %nexcl = pnode(jnode) %nexcl + 1
                    ALLOCATE (current_excl%next,STAT=ios)
                    IF ( ios /= 0 ) CALL stop_memory ( 'exclusion', &
                         'current_excl%next', 0 )
                    current_excl => current_excl%next
                 END IF
              END DO
           END IF DISTIF
           ll_dist => ll_dist%next
        END DO DISTL

! doing the group constraints
        ll_g3x3 => molecule(imol) %ll_g3x3_const
        G3X3L: DO ig3x3 = 1, n3x3_mol
           bend(1) = ll_g3x3%p1%iatom
           bend(2) = ll_g3x3%p2%iatom
           bend(3) = ll_g3x3%p3%iatom

! checking to see if jatom is bonded
           DO ii = 1, 3
              IF (jatom==bend(ii)) THEN
                 match = .TRUE.
                 EXIT
              END IF
           END DO
           G3X3IF: IF (match) THEN
              DO ii = 1, 3
                 CALL checklist(pnode(jnode),match,jatom,bend(ii))
                 IF ( .NOT. match) THEN
                    i = getnode(molecule,bend(ii))
                    current_excl%p => pnode(i) %p
                    pnode(jnode) %nexcl = pnode(jnode) %nexcl + 1
                    ALLOCATE (current_excl%next,STAT=ios)
                    IF ( ios /= 0 ) CALL stop_memory ( 'exclusion', &
                         'current_excl%next', 0 )
                    current_excl => current_excl%next
                 END IF
              END DO
           END IF G3X3IF
           ll_g3x3 => ll_g3x3%next
        END DO G3X3L
        ll_atom => ll_atom%next
     END DO ALOOP
  END DO MLOOP

! deallocating
  DEALLOCATE (bond,STAT=ios)
  IF ( ios /= 0 ) CALL stop_memory ( 'exclusion','bond')
  DEALLOCATE (bend,STAT=ios)
  IF ( ios /= 0 ) CALL stop_memory ( 'exclusion','bend')
  DEALLOCATE (torsion,STAT=ios)
  IF ( ios /= 0 ) CALL stop_memory ( 'exclusion','torsion')

END SUBROUTINE exclusion

!!*****
!******************************************************************************
!!****** linklists/checklist [1.0] *
!!
!!   NAME
!!     checklist
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE checklist(pnode,match,jatom,ipart)

! checks to see if excluded particle is already in the list

  IMPLICIT NONE

! Arguments
  TYPE (particle_node_type ), INTENT ( IN ) :: pnode
  INTEGER, INTENT ( IN ) :: ipart, jatom
  LOGICAL, INTENT ( OUT ) :: match

! Locals
  TYPE ( linklist_exclusion ), POINTER :: ll_exclude
  INTEGER :: i

!------------------------------------------------------------------------------

! initialize match
  match = .FALSE.
!BIGIF:  IF ((jatom==ipart) .OR. (ipart<jatom)) THEN
  BIGIF: IF (jatom==ipart) THEN
     match = .TRUE.
  ELSE
     ll_exclude => pnode%ex
     LOOP: DO i = 1, pnode%nexcl
        IF (ipart==ll_exclude%p%iatom) THEN
           match = .TRUE.
           EXIT LOOP
        END IF
        ll_exclude => ll_exclude%next
     END DO LOOP
  END IF BIGIF

END SUBROUTINE checklist

!!*****
!******************************************************************************

FUNCTION getnode ( molecule, atom ) RESULT ( node )

  IMPLICIT NONE

! Return value
  INTEGER :: node

! Arguments
  TYPE (molecule_structure_type ), DIMENSION ( : ), INTENT ( IN ) :: molecule
  INTEGER, INTENT ( IN ) :: atom

! Locals
  INTEGER :: imol

!------------------------------------------------------------------------------

  node = 0
  DO imol = 1, SIZE ( molecule )
     IF ( molecule ( imol ) % atombase + molecule ( imol ) % natoms_mol &
          > atom ) THEN

        IF ( molecule ( imol ) % atombase > atom ) &
             CALL stop_program ( "getnode", "error in atombase" )

        node = node + ( atom - molecule ( imol ) % atombase + 1 )

        EXIT
     END IF
     node = node + molecule(imol) % natoms_mol
  END DO

END FUNCTION getnode

!******************************************************************************
!!****** linklists/bond_link_list [1.0] *
!!
!!   NAME
!!     bond_link_list
!!
!!   FUNCTION
!!     Generates linked lists for intramolecular dependencies
!!
!!   AUTHOR
!!     CJM & JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (26-01-2002) : index pointer to bond list
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE bond_link_list ( pnode, molecule, mol_setup, bond_parm )

  IMPLICIT NONE

! Arguments
  TYPE ( particle_node_type ), INTENT ( INOUT ), DIMENSION ( : ), TARGET :: &
       pnode
  TYPE ( molecule_structure_type ), INTENT ( INOUT ), DIMENSION ( : ) :: &
       molecule
  TYPE ( molecule_type ), INTENT ( IN ) :: mol_setup ( : )
  TYPE ( bond_parameters_type ), DIMENSION ( : ), POINTER :: bond_parm

! Locals
  INTEGER :: imol, jmol, nmol, ibond, nb, atombase, pTYPE
  INTEGER :: index ( 2 ), ios
  TYPE ( linklist_bonds ), POINTER :: current_ll_bonds

!------------------------------------------------------------------------------

  ALLOCATE (current_ll_bonds,STAT=ios)
  IF ( ios /= 0 ) CALL stop_memory ( 'bond_link_list', 'current_ll_bonds', 0 )

  bonds = .FALSE.
  nmol = size(molecule)
  atombase = 0
  DO imol = 1, nmol
     jmol = molecule(imol) %moltype
     nb = mol_setup(jmol) %molpar%nbonds
     bonds = bonds .OR. (nb /= 0 )
     molecule(imol) %nbonds_mol = nb
     molecule(imol) %ll_bonds => current_ll_bonds
     DO ibond = 1, nb

!
! index(1) and index(2) are the part array indices
! of the first and second atoms in
! the ibond bond pair. pTYPE is the bond type
!
        index(1) = atombase + mol_setup(jmol) %molpar%bonds(1,ibond)
        index(2) = atombase + mol_setup(jmol) %molpar%bonds(2,ibond)
        pTYPE = mol_setup(jmol) %molpar%bonds(3,ibond)
        current_ll_bonds%p1 => pnode(index(1)) %p
        current_ll_bonds%p2 => pnode(index(2)) %p
        current_ll_bonds%index = index
        current_ll_bonds%bond_param => bond_parm(pTYPE)
        ALLOCATE (current_ll_bonds%next,STAT=ios)
        IF ( ios /= 0 ) &
             CALL stop_memory ( 'bond_link_list', 'current_ll_bonds%next', 0 )
        current_ll_bonds => current_ll_bonds%next
     END DO
     atombase = atombase + mol_setup(jmol) %molpar%natom
  END DO

END SUBROUTINE bond_link_list

!!*****
!******************************************************************************
!!****** linklists/bend_link_list [1.0] *
!!
!!   NAME
!!     bend_link_list
!!
!!   FUNCTION
!!     Generates linked lists for intramolecular dependencies
!!
!!   AUTHOR
!!     CJM & JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (26-01-2002) : index pointer to bend list
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE bend_link_list ( pnode, molecule, mol_setup, bend_parm )

  IMPLICIT NONE

! Arguments
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( INOUT ), TARGET :: &
       pnode
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( INOUT ) :: &
       molecule
  TYPE ( molecule_type ), INTENT ( IN ) :: mol_setup ( : )
  TYPE ( bend_parameters_type ), DIMENSION (:), POINTER :: bend_parm

! Locals
  INTEGER :: imol, jmol, nmol, ibend, nb, atombase, pTYPE
  INTEGER :: index ( 3 ), ios
  TYPE ( linklist_bends ), POINTER :: current_ll_bends

!------------------------------------------------------------------------------

  ALLOCATE (current_ll_bends,STAT=ios)
  IF ( ios /= 0 ) CALL stop_memory ( 'bend_link_list', 'current_ll_bends', 0 )

  bends = .FALSE.
  nmol = size(molecule)
  atombase = 0
  DO imol = 1, nmol
     jmol = molecule(imol) %moltype
     nb = mol_setup(jmol) %molpar%nbends
     bends = bends .OR. (nb /= 0 )
     molecule(imol) %nbends_mol = nb
     molecule(imol) %ll_bends => current_ll_bends
     DO ibend = 1, nb

!
! index(1),index(2) and index(3) are the part array indices of the
!    first, second and third atoms in
!    the ibend bend trio. pTYPE is the bend type
!
        index(1) = atombase + mol_setup(jmol) %molpar%bends(1,ibend)
        index(2) = atombase + mol_setup(jmol) %molpar%bends(2,ibend)
        index(3) = atombase + mol_setup(jmol) %molpar%bends(3,ibend)
        pTYPE = mol_setup(jmol) %molpar%bends(4,ibend)
        current_ll_bends%p1 => pnode(index(1)) %p
        current_ll_bends%p2 => pnode(index(2)) %p
        current_ll_bends%p3 => pnode(index(3)) %p
        current_ll_bends%bend_param => bend_parm(pTYPE)
        current_ll_bends%index = index
        ALLOCATE (current_ll_bends%next,STAT=ios)
        IF ( ios /= 0 ) &
             CALL stop_memory ( 'bend_link_list', 'current_ll_bends%next', 0 )
        current_ll_bends => current_ll_bends%next
     END DO
     atombase = atombase + mol_setup(jmol) %molpar%natom
  END DO

END SUBROUTINE bend_link_list

!!*****
!******************************************************************************
!!****** linklists/torsion_link_list [1.0] *
!!
!!   NAME
!!     torsion_link_list
!!
!!   FUNCTION
!!     Generates linked lists for intramolecular dependencies
!!
!!   AUTHOR
!!     CJM & JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (26-01-2002) : index pointer to torsion list
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE torsion_link_list ( pnode, molecule, mol_setup, torsion_parm )

  IMPLICIT NONE

! Arguments
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( INOUT ), TARGET :: &
       pnode
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( INOUT ) :: &
       molecule
  TYPE ( molecule_type ), INTENT ( IN ) :: mol_setup ( : )
  TYPE ( torsion_parameters_type ), DIMENSION (:), POINTER :: torsion_parm

! Locals
  INTEGER :: imol, jmol, nmol, itorsion, nb, atombase, pTYPE
  INTEGER :: index ( 4 ), ios
  TYPE (linklist_torsions), POINTER :: current_ll_torsions

!------------------------------------------------------------------------------

  ALLOCATE (current_ll_torsions,STAT=ios)
  IF ( ios /= 0 ) &
       CALL stop_memory ( 'torsion_link_list', 'current_ll_torsions', 0 )

  torsions = .FALSE.
  nmol = size(molecule)
  atombase = 0
  DO imol = 1, nmol
     jmol = molecule(imol) %moltype
     nb = mol_setup(jmol) %molpar%ntorsions
     torsions = torsions .OR. (nb /= 0 )
     molecule(imol) %ntorsions_mol = nb
     molecule(imol) %ll_torsions => current_ll_torsions
     DO itorsion = 1, nb

!
! index(1),index(2),index(3) and index(4) are the part array indices of the
! first,second,third and fourth atoms in
! the itorsion torsion quartet. pTYPE is the torsion type
!

        index(1) = atombase + mol_setup(jmol) %molpar%torsions(1,itorsion)
        index(2) = atombase + mol_setup(jmol) %molpar%torsions(2,itorsion)
        index(3) = atombase + mol_setup(jmol) %molpar%torsions(3,itorsion)
        index(4) = atombase + mol_setup(jmol) %molpar%torsions(4,itorsion)
        pTYPE = mol_setup(jmol) %molpar%torsions(5,itorsion)
        current_ll_torsions%p1 => pnode(index(1)) %p
        current_ll_torsions%p2 => pnode(index(2)) %p
        current_ll_torsions%p3 => pnode(index(3)) %p
        current_ll_torsions%p4 => pnode(index(4)) %p
        current_ll_torsions%torsion_param => torsion_parm(pTYPE)
        current_ll_torsions%index = index
        ALLOCATE (current_ll_torsions%next,STAT=ios)
        IF ( ios /= 0 ) CALL stop_memory ( 'torsion_link_list', &
             'current_ll_torsions%next', 0 )
        current_ll_torsions => current_ll_torsions%next
     END DO
     atombase = atombase + mol_setup(jmol) %molpar%natom
  END DO

END SUBROUTINE torsion_link_list

!!*****
!******************************************************************************
!!****** linklists/distconst_link_list [1.0] *
!!
!!   NAME
!!     distconst_link_list
!!
!!   FUNCTION
!!     Generates linked lists for intramolecular dependencies
!!
!!   AUTHOR
!!     CJM & JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (26-01-2002) : index pointer to distcon list
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE distconst_link_list ( pnode, molecule, mol_setup, &
     constraint_distance )

  IMPLICIT NONE

! Arguments
  TYPE ( particle_node_type ), INTENT ( INOUT ), DIMENSION ( : ), TARGET :: &
       pnode
  TYPE ( molecule_structure_type ), INTENT ( INOUT ), DIMENSION ( : ) :: &
       molecule
  TYPE ( molecule_type ), INTENT ( IN ) :: mol_setup ( : )
  REAL ( dbl ), DIMENSION ( : ), POINTER :: constraint_distance

! Locals
  INTEGER :: imol, jmol, nmol, idcon, nb, atombase, pTYPE
  INTEGER :: index ( 2 ), ios
  TYPE ( linklist_dist_constraints ), POINTER :: current_ll_dconst

!------------------------------------------------------------------------------

  ALLOCATE (current_ll_dconst,STAT=ios)
  IF ( ios /= 0 ) &
       CALL stop_memory ( 'distconst_link_list', 'current_ll_dconst', 0 )

  dist_constraints = .FALSE.
  nmol = size(molecule)
  atombase = 0
  DO imol = 1, nmol
     jmol = molecule(imol) %moltype
     nb = mol_setup(jmol) %molpar%ndcon
     dist_constraints = dist_constraints .OR. (nb /= 0 )
     molecule(imol) %ndistconst_mol = nb
     molecule(imol) %ll_dist_const => current_ll_dconst
     DO idcon = 1, nb

!
! index(1) and index(2) are the part array indices of the
! first and second atoms in
! the iconst const. pair. pTYPE is the constraint type
!
        index(1) = atombase + mol_setup(jmol) %molpar%ndc(1,idcon)
        index(2) = atombase + mol_setup(jmol) %molpar%ndc(2,idcon)
        pTYPE = mol_setup(jmol) %molpar%ndc(3,idcon)
        current_ll_dconst%p1 => pnode(index(1)) %p
        current_ll_dconst%p2 => pnode(index(2)) %p
        current_ll_dconst%distance => constraint_distance(pTYPE)
        current_ll_dconst%lambda = 0._dbl
        current_ll_dconst%fc ( : ) = 0._dbl
        ALLOCATE (current_ll_dconst%next,STAT=ios)
        IF ( ios /= 0 ) CALL stop_memory ( 'distconst_link_list', &
             'current_ll_dconst%next', 0 )

        current_ll_dconst => current_ll_dconst % next
     END DO

     atombase = atombase + mol_setup ( jmol ) % molpar % natom
  END DO

END SUBROUTINE distconst_link_list

!!*****
!******************************************************************************
!!****** linklists/g3x3_link_list [1.0] *
!!
!!   NAME
!!     g3x3_link_list
!!
!!   FUNCTION
!!     Generates linked lists for intramolecular dependencies
!!
!!   AUTHOR
!!     CJM & JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (26-01-2002) : index pointer to distcon list
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE g3x3_link_list ( pnode, molecule, mol_setup, constraint_distance )

  IMPLICIT NONE

! Arguments
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( INOUT ), TARGET :: &
       pnode
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( INOUT ) :: &
       molecule
  TYPE ( molecule_type ), INTENT ( IN ) :: mol_setup ( : )
  REAL ( dbl ), DIMENSION ( : ), POINTER :: constraint_distance
! Locals
  INTEGER :: imol, jmol, nmol, ig3x3, nb, atombase, pTYPE ( 3 )
  INTEGER :: index ( 3 ), ios
  TYPE ( linklist_g3x3_constraints ), POINTER :: current_ll_g3x3

!------------------------------------------------------------------------------

  ALLOCATE ( current_ll_g3x3, STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'g3x3_link_list', 'current_ll_g3x3', 0 )

  g3x3_constraints = .FALSE.
  nmol = size ( molecule )
  atombase = 0

  DO imol = 1, nmol
     jmol = molecule(imol) % moltype
     nb = mol_setup(jmol) % molpar % n3x3con
     g3x3_constraints = g3x3_constraints .OR. ( nb /= 0 )
     molecule(imol) %ng3x3const_mol = nb
     molecule(imol) %ll_g3x3_const => current_ll_g3x3
     DO ig3x3 = 1, nb

!
!    index(1),index(2) and index(3) are the part array indices of the
!    first,second and third atoms in
!    the ibend bend trio. pTYPE(:) are the constraint types
!
        index(1) = atombase + mol_setup(jmol) %molpar%n3x3c(1,ig3x3)
        index(2) = atombase + mol_setup(jmol) %molpar%n3x3c(2,ig3x3)
        index(3) = atombase + mol_setup(jmol) %molpar%n3x3c(3,ig3x3)
        pTYPE ( 1) = mol_setup(jmol) %molpar%n3x3c(4,ig3x3)
        pTYPE ( 2) = mol_setup(jmol) %molpar%n3x3c(5,ig3x3)
        pTYPE ( 3) = mol_setup(jmol) %molpar%n3x3c(6,ig3x3)
        current_ll_g3x3%p1 => pnode(index(1)) % p
        current_ll_g3x3%p2 => pnode(index(2)) % p
        current_ll_g3x3%p3 => pnode(index(3)) % p
        current_ll_g3x3%d12 => constraint_distance(pTYPE ( 1))
        current_ll_g3x3%d13 => constraint_distance(pTYPE ( 2))
        current_ll_g3x3%d23 => constraint_distance(pTYPE ( 3))
        current_ll_g3x3%lambda ( : ) = 0.0_dbl
        current_ll_g3x3%fc1 ( : ) = 0.0_dbl
        current_ll_g3x3%fc2 ( : ) = 0.0_dbl
        current_ll_g3x3%fc3 ( : ) = 0.0_dbl

        ALLOCATE ( current_ll_g3x3 % next, STAT = ios )
        IF ( ios /= 0 ) &
             CALL stop_memory ( 'g3x3_link_list', 'current_ll_g3x3%next', 0 )

        current_ll_g3x3 => current_ll_g3x3 % next
     END DO

     atombase = atombase + mol_setup ( jmol ) % molpar % natom
  END DO

END SUBROUTINE g3x3_link_list

!!*****
!******************************************************************************
!!****** linklists/g4x6_link_list [1.0] *
!!
!!   NAME
!!     g4x6_link_list
!!
!!   FUNCTION
!!     Generates linked lists for intramolecular dependencies
!!
!!   AUTHOR
!!     CJM & JGH
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE g4x6_link_list ( pnode, molecule, mol_setup, constraint_distance )

  IMPLICIT NONE

! Arguments
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( INOUT ), TARGET :: &
       pnode
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( INOUT ) :: &
       molecule
  TYPE ( molecule_type ), INTENT ( IN ) :: mol_setup ( : )
  REAL ( dbl ), DIMENSION ( : ), POINTER :: constraint_distance
! Locals
  INTEGER :: imol, jmol, nmol, ig4x6, nb, atombase, pTYPE ( 6 )
  INTEGER :: index ( 4 ), ios
  TYPE ( linklist_g4x6_constraints ), POINTER :: current_ll_g4x6

!------------------------------------------------------------------------------

  ALLOCATE ( current_ll_g4x6, STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'g4x6_link_list', 'current_ll_g4x6', 0 )

  g4x6_constraints = .FALSE.
  nmol = size ( molecule )
  atombase = 0

  DO imol = 1, nmol
     jmol = molecule(imol) % moltype
     nb = mol_setup(jmol) % molpar % n3x3con
     g4x6_constraints = g4x6_constraints .OR. ( nb /= 0 )
     molecule(imol) %ng4x6const_mol = nb
     molecule(imol) %ll_g4x6_const => current_ll_g4x6
     DO ig4x6 = 1, nb

!
!    index(1),index(2),index(3) and index(4) are the part array indices of the
!    first,second,third and fourth atoms
!    pTYPE(:) are the constraint types
!
        index(1) = atombase + mol_setup(jmol) %molpar%n4x6c(1,ig4x6)
        index(2) = atombase + mol_setup(jmol) %molpar%n4x6c(2,ig4x6)
        index(3) = atombase + mol_setup(jmol) %molpar%n4x6c(3,ig4x6)
        index(4) = atombase + mol_setup(jmol) %molpar%n4x6c(4,ig4x6)
        pTYPE ( 1) = mol_setup(jmol) %molpar%n4x6c(5,ig4x6)
        pTYPE ( 2) = mol_setup(jmol) %molpar%n4x6c(6,ig4x6)
        pTYPE ( 3) = mol_setup(jmol) %molpar%n4x6c(7,ig4x6)
        pTYPE ( 4) = mol_setup(jmol) %molpar%n4x6c(8,ig4x6)
        pTYPE ( 5) = mol_setup(jmol) %molpar%n4x6c(9,ig4x6)
        pTYPE ( 6) = mol_setup(jmol) %molpar%n4x6c(10,ig4x6)
        current_ll_g4x6%p1 => pnode(index(1)) % p
        current_ll_g4x6%p2 => pnode(index(2)) % p
        current_ll_g4x6%p3 => pnode(index(3)) % p
        current_ll_g4x6%p4 => pnode(index(4)) % p
        current_ll_g4x6%d12 => constraint_distance(pTYPE ( 1))
        current_ll_g4x6%d13 => constraint_distance(pTYPE ( 2))
        current_ll_g4x6%d14 => constraint_distance(pTYPE ( 3))
        current_ll_g4x6%d23 => constraint_distance(pTYPE ( 4))
        current_ll_g4x6%d24 => constraint_distance(pTYPE ( 5))
        current_ll_g4x6%d34 => constraint_distance(pTYPE ( 6))
        current_ll_g4x6%lambda ( : ) = 0.0_dbl
        current_ll_g4x6%fc1 ( : ) = 0.0_dbl
        current_ll_g4x6%fc2 ( : ) = 0.0_dbl
        current_ll_g4x6%fc3 ( : ) = 0.0_dbl
        current_ll_g4x6%fc4 ( : ) = 0.0_dbl
        current_ll_g4x6%fc5 ( : ) = 0.0_dbl
        current_ll_g4x6%fc6 ( : ) = 0.0_dbl

        ALLOCATE ( current_ll_g4x6 % next, STAT = ios )
        IF ( ios /= 0 ) &
             CALL stop_memory ( 'g4x6_link_list', 'current_ll_g4x6%next', 0 )

        current_ll_g4x6 => current_ll_g4x6 % next
     END DO

     atombase = atombase + mol_setup ( jmol ) % molpar % natom
  END DO

END SUBROUTINE g4x6_link_list

!!*****
!******************************************************************************

END MODULE linklists

!******************************************************************************
