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
!!     CJM & JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE dft_types
  
  USE kinds, ONLY : dbl
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: dft_control_type, wave_control_type
  
  !! dft_control_type = general DFT
  !! wave_control_type = options/control of the DFT/PW/PP part
  
  TYPE wave_control_type
     REAL ( dbl ) :: e_cutoff_wf, e_cutoff_dual
     INTEGER :: xc_mesh ( 3 )
     
     INTEGER :: niter_elc
     INTEGER :: wf_print_frequency
     CHARACTER ( LEN = 80 ) :: wf_optimisation_method
     REAL ( dbl ) :: wf_timestep
     REAL ( dbl ) :: wf_emass
     LOGICAL :: wf_subspacerotation_tperform
     INTEGER :: wf_subspacerotation_frequency
     INTEGER :: wf_subspacerotation_first
  END TYPE wave_control_type
  
  !! niter_elc = maximum number of electronic iterations before ionic move
  TYPE dft_control_type
     CHARACTER ( LEN = 40 ) :: functional
     REAL ( dbl ) :: density_cut
     REAL ( dbl ) :: gradient_cut
     REAL ( dbl ) :: nelectrons
     REAL ( dbl ) :: fixedspinmoment_amount
     LOGICAL :: broadening_tperform, fixedspinmoment
     CHARACTER ( LEN = 40 ) :: broadening_method
     REAL ( dbl ) :: broadening_width
     
     INTEGER :: nstates, nspins
     REAL ( dbl ) :: charge
     
     TYPE ( wave_control_type ) :: wave_control
  END TYPE dft_control_type
  
END MODULE dft_types

!!*****
!******************************************************************************
