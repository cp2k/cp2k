!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/global_types [1.0] *
!!
!!   NAME
!!     global_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     - print keys, basis_set_file name and potential_file_name added to the
!!       global type (27.02.2001, MK)
!!
!!   SOURCE
!******************************************************************************

MODULE global_types

  USE kinds, ONLY : dbl

  USE print_keys, ONLY: print_key_type

  IMPLICIT NONE

  PRIVATE

  TYPE global_environment_type
     TYPE(print_key_type), POINTER :: print
     CHARACTER ( LEN = 200 ) :: basis_set_file_name,&
                                input_file_name,&
                                potential_file_name,&
                                restart_file_name
     CHARACTER ( LEN = 200 ) :: program_name
     CHARACTER ( LEN = 200 ) :: project_name
     CHARACTER ( LEN = 200 ) :: default_fft_library
     LOGICAL :: ionode        !! this processor performs global output
     INTEGER :: group         !! group handle for all processors
     INTEGER :: source        !! processor id for global output node
     INTEGER :: num_pe        !! total number of processors
     INTEGER :: mepos         !! position number for local processor
     INTEGER :: print_level   !! how much information is written
     INTEGER :: idum          !! random number seed
     INTEGER :: scr           !! output unit
  END TYPE global_environment_type

  PUBLIC :: global_environment_type

END MODULE global_types

!!*****
!******************************************************************************
