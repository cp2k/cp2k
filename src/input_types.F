!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE input_types
  
  IMPLICIT NONE
  
  TYPE setup_parameters_type
     INTEGER                :: nmol, natoms, natom_type
     INTEGER                :: perd ( 3 )
     LOGICAL                :: symmetry
     CHARACTER ( LEN = 20 ) :: set_file_name
     CHARACTER ( LEN = 20 ) :: input_file_name
     CHARACTER ( LEN = 10 ) :: run_type
     CHARACTER ( LEN = 10 ) :: unit_type
  END TYPE setup_parameters_type
  
  PRIVATE
  PUBLIC :: setup_parameters_type
  
END MODULE input_types
