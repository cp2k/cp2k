!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/dft_types [1.0] *
!!
!!   NAME
!!     dft_types
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

MODULE dft_types
  
  USE band, ONLY : band_structure_type
  USE kinds, ONLY : dbl
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: dft_control_type, wave_control_type, qs_control_type
  
  TYPE qs_control_type
     CHARACTER ( LEN = 10 ) :: qs_method
  END TYPE qs_control_type

  TYPE wave_control_type
     REAL ( dbl ) :: e_cutoff_wf
     REAL ( dbl ) :: e_cutoff_dual
  END TYPE wave_control_type
  
  TYPE dft_control_type
     CHARACTER ( LEN = 40 ) :: functional
     REAL ( dbl ) :: density_cut
     REAL ( dbl ) :: gradient_cut
     INTEGER :: xc_mesh ( 3 )
     INTEGER :: nspins
     
     TYPE ( wave_control_type ) :: wave_control
     TYPE ( qs_control_type ) :: qs_control

     TYPE ( band_structure_type ) :: band
  END TYPE dft_control_type
  
END MODULE dft_types

!!*****
!******************************************************************************
