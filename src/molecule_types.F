!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/molecule_types [1.0] *
!!
!!   NAME
!!     molecule_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM & JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (26-01-2002) : Store bonds, bends etc in tables not in matrices
!!
!!   SOURCE
!******************************************************************************

MODULE molecule_types

  USE kinds, ONLY : dbl
  USE particle_types, ONLY : particle_type
  USE qs_neighbor_list_types, ONLY: neighbor_list_type

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: molecule_structure_type, linklist_atoms, bond_parameters_type, &
       linklist_bonds, bend_parameters_type, linklist_bends,  &
       torsion_parameters_type, linklist_torsions, linklist_dist_constraints, &
       linklist_g3x3_constraints, linklist_g4x6_constraints, linklist_exclusion, &
       particle_node_type, localization_list
  PUBLIC :: particle_type
  PUBLIC :: molecule_type, molecule_parameter_type
  PUBLIC :: intra_parameters_type
  PUBLIC :: topology_type

  ! public procedures
  PUBLIC :: init_molecule_type

  TYPE linklist_exclusion
     TYPE ( particle_type ), POINTER :: p
     TYPE ( linklist_exclusion ), POINTER :: next
  END TYPE linklist_exclusion

  TYPE localization_list
     LOGICAL :: all
     INTEGER :: nloc
     INTEGER, DIMENSION ( : ), POINTER :: atomlist
  END TYPE localization_list

  TYPE particle_node_type
     TYPE ( particle_type ), POINTER :: p
     INTEGER :: in_mol, iat_in_mol
     INTEGER :: nneighbor, nexcl
     TYPE ( neighbor_list_type ), POINTER :: nl
     TYPE ( linklist_exclusion ), POINTER :: ex
     TYPE ( localization_list ), POINTER :: llist
  END TYPE particle_node_type

  TYPE linklist_atoms
     TYPE ( particle_node_type ), POINTER :: part
     TYPE ( linklist_atoms ), POINTER :: next
  END TYPE linklist_atoms

  TYPE bond_parameters_type
     CHARACTER ( LEN = 20 ) :: type
     REAL ( dbl ) :: r0, k
  END TYPE bond_parameters_type

  TYPE linklist_bonds
     TYPE ( particle_type ), POINTER :: p1, p2
     TYPE ( bond_parameters_type ), POINTER :: bond_param
     TYPE ( linklist_bonds), POINTER :: next
     INTEGER, DIMENSION ( 2 ) :: index
  END TYPE linklist_bonds

  TYPE bend_parameters_type
     CHARACTER ( LEN = 20 ) :: type
     REAL ( dbl ) :: theta0, k
  END TYPE bend_parameters_type

  TYPE linklist_bends
     TYPE ( particle_type ), POINTER :: p1, p2, p3
     TYPE ( bend_parameters_type ), POINTER :: bend_param
     TYPE ( linklist_bends ), POINTER :: next
     INTEGER, DIMENSION ( 3 ) :: index
  END TYPE linklist_bends

  TYPE torsion_parameters_type
     CHARACTER ( LEN = 20 ) :: type
     REAL ( dbl ), DIMENSION (0:3) :: a
  END TYPE torsion_parameters_type

  TYPE linklist_torsions
     TYPE ( particle_type ), POINTER :: p1, p2, p3, p4
     TYPE ( torsion_parameters_type ), POINTER :: torsion_param
     TYPE ( linklist_torsions ), POINTER :: next
     INTEGER, DIMENSION ( 4 ) :: index
  END TYPE linklist_torsions

  TYPE linklist_dist_constraints
     TYPE ( particle_type ), POINTER :: p1, p2
     REAL ( dbl ), POINTER :: distance
     REAL ( dbl ) :: del_lam, lambda, sigma
     REAL ( dbl ), DIMENSION (3) :: rold1, rold2, fc
     REAL ( dbl ) :: scale, scale_old
     TYPE ( linklist_dist_constraints ), POINTER :: next
  END TYPE linklist_dist_constraints

  TYPE linklist_g3x3_constraints
     TYPE ( particle_type ), POINTER :: p1, p2, p3
     REAL ( dbl ), POINTER :: d12, d13, d23
     REAL ( dbl ), DIMENSION( 3 ) :: fc1, fc2, fc3
     REAL ( dbl ), DIMENSION( 3 ) :: lambda, del_lam
     REAL ( dbl ), DIMENSION( 3 ) :: rold1, rold2, rold3
     REAL ( dbl ), DIMENSION( 3 ) :: vold1, vold2, vold3
     REAL ( dbl ), DIMENSION( 3 ) :: v1, v2, v3
     REAL ( dbl ) :: scale, scale_old
     TYPE ( linklist_g3x3_constraints ), POINTER :: next
  END TYPE linklist_g3x3_constraints

  TYPE linklist_g4x6_constraints
     TYPE ( particle_type ), POINTER :: p1, p2, p3, p4
     REAL ( dbl ), POINTER :: d12, d13, d14, d23, d24, d34
     REAL ( dbl ), DIMENSION( 3 ) :: fc1, fc2, fc3, fc4, fc5, fc6
     REAL ( dbl ), DIMENSION( 6 ) :: lambda, del_lam
     REAL ( dbl ), DIMENSION( 3 ) :: rold1, rold2, rold3, rold4
     REAL ( dbl ), DIMENSION( 3 ) :: vold1, vold2, vold3, vold4
     REAL ( dbl ), DIMENSION( 3 ) :: v1, v2, v3, v4
     REAL ( dbl ) :: scale, scale_old
     TYPE ( linklist_g4x6_constraints ), POINTER :: next
  END TYPE linklist_g4x6_constraints

  TYPE intra_parameters_type
     TYPE ( bond_parameters_type ), POINTER :: bond_param ( : )
     TYPE ( bend_parameters_type ), POINTER :: bend_param ( : )
     TYPE ( torsion_parameters_type ), POINTER :: torsion_param ( : )
     REAL ( dbl ), POINTER :: constraint_distance ( : )
  END TYPE intra_parameters_type

  TYPE molecule_parameter_type
     CHARACTER ( LEN = 20 ) :: name
     INTEGER :: natom
     INTEGER :: natom_type
     CHARACTER ( LEN = 20 ), POINTER :: aname ( : )
     CHARACTER ( LEN = 20 ), POINTER :: alabel ( : )
     REAL ( dbl ), POINTER :: aweight ( : )
     REAL ( dbl ), POINTER :: acharge ( : )
     INTEGER :: nbonds
     INTEGER, POINTER :: bonds ( :, : )
     INTEGER :: nbends
     INTEGER, POINTER :: bends ( :, : )
     INTEGER :: ntorsions
     INTEGER, POINTER :: torsions ( :, : )
     INTEGER :: ndcon
     INTEGER, POINTER :: ndc ( :, : )
     INTEGER :: n3x3con
     INTEGER, POINTER :: n3x3c ( :, : )
     INTEGER :: n4x6con
     INTEGER, POINTER :: n4x6c ( :, : )
  END TYPE molecule_parameter_type

  TYPE molecule_type
     CHARACTER ( LEN = 20 ) :: file_name
     INTEGER :: num_mol
     CHARACTER ( LEN = 20 ) :: nhcopt
     TYPE ( molecule_parameter_type ) :: molpar
  END TYPE molecule_type

  TYPE topology_type
     INTEGER :: nbond
     INTEGER, DIMENSION (:,:), POINTER :: bond
     INTEGER :: nbend
     INTEGER, DIMENSION (:,:), POINTER :: bend
     INTEGER :: ntors
     INTEGER, DIMENSION (:,:), POINTER :: torsion
  END TYPE topology_type

  TYPE molecule_structure_type
     CHARACTER ( LEN = 30 ) :: mol_name ! name of molecule
     INTEGER :: nstates ! number of electronic states 
     INTEGER, POINTER :: states ( : ) ! indicies of Kohn-Sham state for molecule
     INTEGER :: moltype ! pointer to molecular type
     INTEGER :: imol ! position counter
     INTEGER :: atombase ! position counter of first atom in molecule
     INTEGER :: natoms_mol ! # of atoms in molecule
     INTEGER :: nvt_num ! thermostat number that it is connected.
! Used for constraints
     TYPE (linklist_atoms), POINTER :: ll_atoms
     INTEGER :: nbonds_mol ! # of bonds in molecule
     TYPE (linklist_bonds), POINTER :: ll_bonds
     INTEGER :: nbends_mol ! # of bends in molecule
     TYPE (linklist_bends), POINTER :: ll_bends
     INTEGER :: ntorsions_mol ! # of torsions in molecule
     TYPE (linklist_torsions), POINTER :: ll_torsions
     INTEGER :: ndistconst_mol ! # of distance constraints in molecule
     TYPE (linklist_dist_constraints), POINTER :: ll_dist_const
     INTEGER :: ng3x3const_mol ! # of 3x3 group constraints in molecule
     TYPE (linklist_g3x3_constraints), POINTER :: ll_g3x3_const
     INTEGER :: ng4x6const_mol ! # of 4x6 group constraints in molecule
     TYPE (linklist_g4x6_constraints), POINTER :: ll_g4x6_const
  END TYPE molecule_structure_type

CONTAINS

  SUBROUTINE init_molecule_parameter_type(molpar)
    TYPE(molecule_parameter_type) :: molpar
    molpar%name="#"
    molpar%natom=0
    molpar%natom_type=0
    molpar%nbonds=0
    molpar%nbends=0
    molpar%ntorsions=0
    molpar%ndcon=0
    molpar%n3x3con=0
    molpar%n4x6con=0
    NULLIFY(molpar%aname)
    NULLIFY(molpar%alabel)
    NULLIFY(molpar%aweight)
    NULLIFY(molpar%acharge) 
    NULLIFY(molpar%bonds)
    NULLIFY(molpar%bends)
    NULLIFY(molpar%torsions)
    NULLIFY(molpar%ndc)
    NULLIFY(molpar%n3x3c)
    NULLIFY(molpar%n4x6c)
  END SUBROUTINE

  SUBROUTINE init_molecule_type(mol)
    TYPE(molecule_type) :: mol
    mol%file_name=""
    mol%nhcopt=""
    mol%num_mol=0
    CALL init_molecule_parameter_type(mol%molpar)
  END SUBROUTINE

END MODULE molecule_types

!******************************************************************************
!!*****
