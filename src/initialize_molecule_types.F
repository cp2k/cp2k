!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/initialize_molecule_types [1.0] *
!!
!!   NAME
!!     initialize_molecule_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE initialize_molecule_types

  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE linklists, ONLY : atom_link_list, exclusion, bond_link_list, &
       bend_link_list, torsion_link_list, distconst_link_list, g3x3_link_list
  USE md, ONLY : simulation_parameters_type
  USE molecule_types, ONLY : molecule_type, intra_parameters_type, &
       particle_node_type, molecule_structure_type
  USE particle_types, ONLY : particle_type, particle_prop_type
  USE util, ONLY : get_limit

  IMPLICIT NONE

  PUBLIC :: initialize_molecule_type
  PRIVATE

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** initialize_molecule_types/initialize_molecule_type [1.0] *
!!
!!   NAME
!!     initialize_molecule_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE initialize_molecule_type ( mol_setup, intra_param, &
     pnode, part, molecule, globenv )

  IMPLICIT NONE
  TYPE ( molecule_type ), DIMENSION ( : ), INTENT ( IN ) :: mol_setup
  TYPE ( intra_parameters_type ), INTENT ( IN ) :: intra_param
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( INOUT ) :: pnode
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( INOUT ), TARGET :: part
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( INOUT ) :: &
       molecule
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Locals
  INTEGER :: imol, i, j, k, atombase, iat, nlim ( 2 ), iw, nmol_type

!------------------------------------------------------------------------------

  iw = globenv % scr
  IF ( globenv % num_pe == 1 ) THEN
     imol = 0
     atombase = 1
     nmol_type = SIZE ( mol_setup )
     DO i = 1, nmol_type
        DO k = 1, mol_setup ( i ) %num_mol
           imol = imol + 1
           molecule(imol) % mol_name = mol_setup(i) % molpar% name
           molecule(imol) % imol = imol
           molecule(imol) % atombase = atombase
           molecule(imol) % moltype = i
           molecule(imol) % natoms_mol = mol_setup(i) % molpar% natom
           molecule(imol) % nbonds_mol = mol_setup(i) % molpar% nbonds
           molecule(imol) % nbends_mol = mol_setup(i) % molpar% nbends
           molecule(imol) % ntorsions_mol = mol_setup(i) % molpar % ntorsions
           molecule(imol) % ndistconst_mol = mol_setup(i) % molpar % ndcon
           molecule(imol) % ng3x3const_mol = mol_setup(i) % molpar % n3x3con
           IF (globenv%ionode .AND. globenv%print_level>9) WRITE (iw, &
                '(A,I10,T40,A,I5,A,T61,A20 )' ) ' Molecule:', imol, &
                ' Type:', i, '     Name:', &
                ADJUSTR ( mol_setup(i) % molpar % name)
           atombase = atombase + mol_setup(i) % molpar % natom
        END DO
     END DO
     DO iat = 1, size(pnode)
        pnode(iat) %p => part(iat)
     END DO

  ELSE

     imol = 0
     atombase = 1
     nmol_type = SIZE ( mol_setup )
     DO i = 1, nmol_type
        nlim = get_limit ( mol_setup ( i ) % num_mol, globenv % num_pe, &
             globenv % mepos )
        j = SUM( mol_setup(1:i) % num_mol )
        DO k = nlim(1), nlim(2)
           imol = imol + 1
           molecule(imol) % mol_name = mol_setup(i) % molpar%name
           molecule(imol) % imol = j + k
           molecule(imol) % atombase = atombase + &
                (k-1)*mol_setup(i) % molpar % natom
           molecule(imol) % moltype = i
           molecule(imol) % natoms_mol = mol_setup(i) % molpar%natom
           molecule(imol) % nbonds_mol = mol_setup(i) % molpar%nbonds
           molecule(imol) % nbends_mol = mol_setup(i) % molpar%nbends
           molecule(imol) % ntorsions_mol = mol_setup(i) % molpar%ntorsions
           molecule(imol) % ndistconst_mol = mol_setup(i) % molpar%ndcon
           molecule(imol) % ng3x3const_mol = mol_setup(i) % molpar%n3x3con
           IF (globenv%ionode .AND. globenv%print_level>9) WRITE (iw, &
                '(A,I10,T40,A,I5,A,T61,A20 )' ) ' Molecule:', imol, &
                ' Type:', i, '     Name:', ADJUSTR ( mol_setup(i) % molpar%name)
        END DO
        atombase = atombase &
             + mol_setup(i) % num_mol * mol_setup(i) % molpar % natom
     END DO

     k = 0
     DO imol = 1, size(molecule)
        atombase = molecule(imol) % atombase - 1
        DO j = 1, molecule(imol) % natoms_mol
           k = k + 1
           pnode ( k ) % p => part ( atombase + j )
        END DO
     END DO
  END IF

!..atom link list
  CALL atom_link_list(molecule,pnode)

!..initialize bonds
  CALL bond_link_list(pnode,molecule,mol_setup,intra_param % bond_param)

!..initialize bends
  CALL bend_link_list(pnode,molecule,mol_setup,intra_param % bend_param)

!..initialize torsions
  CALL torsion_link_list(pnode,molecule,mol_setup,intra_param % torsion_param)

!..initialize distance constraints
  CALL distconst_link_list(pnode,molecule,mol_setup, &
       intra_param % constraint_distance)

!..initialize g3x3 constraints
  CALL g3x3_link_list(pnode,molecule,mol_setup, &
       intra_param % constraint_distance)

!..create exclusion list
  CALL exclusion(molecule,pnode)

END SUBROUTINE initialize_molecule_type

!!*****
!******************************************************************************

END MODULE initialize_molecule_types

!******************************************************************************
