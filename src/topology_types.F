!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/topology_types [1.0] *
!!
!!   NAME
!!     topology_types
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

MODULE topology_types

  USE kinds,                           ONLY: dbl,&
                                             default_string_length
  USE simulation_cell,                 ONLY: cell_type

  IMPLICIT NONE

  TYPE topology_parameters_type
     TYPE ( cell_type ), POINTER               :: cell
     INTEGER                                   :: nmol, natoms, natom_type
     LOGICAL                                   :: coordinate
     CHARACTER (LEN=default_string_length)     :: coord_type
     CHARACTER (LEN=default_string_length)     :: coord_file_name
     LOGICAL                                   :: connectivity
     CHARACTER (LEN=default_string_length)     :: conn_type
     CHARACTER (LEN=default_string_length)     :: conn_file_name
     CHARACTER (LEN=default_string_length),&
               DIMENSION ( : ), POINTER        :: atom_names


  CHARACTER (LEN=20), POINTER :: label_resname(:),label_atmname(:)
  REAL(dbl), DIMENSION(:,:), POINTER :: r
  INTEGER, POINTER :: map_mol_typ(:)
  INTEGER, POINTER :: map_mol_num(:)
  INTEGER, POINTER :: bond_a(:),bond_b(:)
  INTEGER, POINTER :: theta_a(:),theta_b(:),theta_c(:)
  INTEGER, POINTER :: phi_a(:),phi_b(:),phi_c(:),phi_d(:)
  REAL (dbl), POINTER :: atm_charge(:)
  REAL (dbl), POINTER :: atm_mass(:)

  END TYPE topology_parameters_type

  PRIVATE
  PUBLIC :: topology_parameters_type

END MODULE topology_types

!!*****
!******************************************************************************
