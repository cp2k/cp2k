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

  USE kinds,                           ONLY: dp,&
                                             default_string_length

  IMPLICIT NONE

  TYPE setup_parameters_type
     LOGICAL                                   :: symmetry
     CHARACTER ( LEN = default_string_length ) :: run_type
     CHARACTER ( LEN = default_string_length ) :: set_file_name
     CHARACTER ( LEN = default_string_length ) :: unit_type
     LOGICAL                                   :: path_integrals
     REAL (KIND=dp)                              :: charge
     INTEGER                                   :: states
     INTEGER                                   :: nspins
     INTEGER                                   :: multiplicity
     REAL (KIND=dp)                              :: spinmoment
     REAL (KIND=dp)                              :: etemp
     CHARACTER ( LEN = default_string_length ) :: broadening_method
  END TYPE setup_parameters_type

  PRIVATE
  PUBLIC :: setup_parameters_type

END MODULE input_types

!!*****
!******************************************************************************
