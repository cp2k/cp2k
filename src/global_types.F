!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/global_types [1.0] *
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
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE global_types
  
  USE kinds, ONLY : dbl
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: global_environment_type
  
  TYPE global_environment_type
     CHARACTER ( LEN = 200 ) :: input_file_name
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
  
END MODULE global_types

!!*****
!******************************************************************************
