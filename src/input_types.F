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
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE input_types

  IMPLICIT NONE

  TYPE setup_parameters_type
     INTEGER                :: nmol, natoms, natom_type
     INTEGER                :: perd ( 3 )
     LOGICAL                :: symmetry
     CHARACTER ( LEN = 40 ) :: set_file_name
     CHARACTER ( LEN = 40 ) :: input_file_name
     CHARACTER ( LEN = 20 ) :: run_type
     CHARACTER ( LEN = 20 ) :: unit_type
  END TYPE setup_parameters_type

  PRIVATE
  PUBLIC :: setup_parameters_type

END MODULE input_types

!!*****
!******************************************************************************
