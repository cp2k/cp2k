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

  USE kinds, ONLY : dbl

  IMPLICIT NONE

  TYPE setup_parameters_type
     INTEGER                :: nmol, natoms, natom_type
     INTEGER                :: perd ( 3 )
     LOGICAL                :: symmetry
     CHARACTER ( LEN = 40 ) :: set_file_name
     CHARACTER ( LEN = 40 ) :: coord_file_name
     CHARACTER ( LEN = 20 ) :: run_type
     CHARACTER ( LEN = 20 ) :: unit_type
     LOGICAL                :: path_integrals
     REAL ( dbl )           :: charge
     INTEGER                :: states
     INTEGER                :: nspins
     INTEGER                :: multiplicity
     REAL ( dbl )           :: spinmoment
     REAL ( dbl )           :: etemp
     CHARACTER ( LEN = 40 ) :: broadening_method
  END TYPE setup_parameters_type

  PRIVATE
  PUBLIC :: setup_parameters_type

END MODULE input_types

!!*****
!******************************************************************************
