!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/input_types [1.0] *
!!
!!   NAME
!!     input_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM & JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (30.11.2001) : new entries in setup_parameters_type
!!                        change name from input_file_name to coord_...
!!                        added topology file
!!                        added atom_names
!!
!!   SOURCE
!******************************************************************************

MODULE input_types

  USE kinds, ONLY : dbl
  USE simulation_cell,                 ONLY: cell_type

  IMPLICIT NONE

  INTEGER, PARAMETER :: ATOMNAMESLENGTH = 20

  TYPE cell_setup_type
     REAL ( dbl ), DIMENSION ( 3, 3 ) :: hmat
     CHARACTER (LEN=20) :: poisson_solver
     INTEGER, DIMENSION (3) :: perd
  END TYPE cell_setup_type

  TYPE setup_parameters_type
     INTEGER                :: nmol, natoms, natom_type
     LOGICAL                :: symmetry
     CHARACTER ( LEN = 40 ) :: set_file_name
     LOGICAL                :: coordinate
     CHARACTER ( LEN = 40 ) :: coordinate_type
     CHARACTER ( LEN = 40 ) :: coord_file_name
     LOGICAL                :: topology
     CHARACTER ( LEN = 40 ) :: topology_type
     CHARACTER ( LEN = 40 ) :: topology_file_name
     CHARACTER ( LEN = 20 ) :: unit_type
     LOGICAL                :: path_integrals
     REAL ( dbl )           :: charge
     INTEGER                :: states
     INTEGER                :: nspins
     INTEGER                :: multiplicity
     REAL ( dbl )           :: spinmoment
     REAL ( dbl )           :: etemp
     CHARACTER ( LEN = 40 ) :: broadening_method
     CHARACTER ( LEN = ATOMNAMESLENGTH ), DIMENSION ( : ), POINTER :: atom_names
     TYPE ( cell_setup_type ) :: simulation_cell
     TYPE ( cell_setup_type ) :: reference_cell
     TYPE ( cell_type ), POINTER :: cell
  END TYPE setup_parameters_type

  PRIVATE
  PUBLIC :: setup_parameters_type, cell_setup_type, ATOMNAMESLENGTH

END MODULE input_types

!!*****
!******************************************************************************
