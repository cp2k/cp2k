!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE fit_parameters_types

  USE kinds,                           ONLY: dbl

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: linklist_real_type, linklist_parameters_type, &
            fit_parameter_type, fit_info_type, fit_averages_type

  TYPE linklist_real_type
     REAL ( dbl ), POINTER :: real_value
     TYPE ( linklist_real_type ), POINTER :: next
  END TYPE linklist_real_type

  TYPE linklist_parameters_type
     TYPE ( linklist_real_type ), POINTER :: lleps
     TYPE ( linklist_real_type ), POINTER :: llsig
     TYPE ( linklist_real_type ), POINTER :: lla
     TYPE ( linklist_real_type ), POINTER :: llb
     TYPE ( linklist_real_type ), POINTER :: llc
  END TYPE linklist_parameters_type
 
  TYPE fit_parameter_type
     CHARACTER ( len = 10 ) :: datatype, outunit, restart
     CHARACTER ( len = 10 ) :: fit_type, pot_type
     INTEGER :: max_iter, nsteps, iskip, iwrite, idump
     REAL ( dbl ) :: max_tol
  END TYPE fit_parameter_type
 
  TYPE fit_info_type
     REAL ( dbl ), DIMENSION ( 3 ) :: f_ref, f_fit, r
  END TYPE fit_info_type

  TYPE fit_averages_type
    REAL (dbl), POINTER, DIMENSION(:,:) :: ave_sig, ave_eps
    REAL (dbl), POINTER, DIMENSION(:,:) :: ave_a, ave_b, ave_c
  END TYPE fit_averages_type
                                                                                          

END MODULE fit_parameters_types

!******************************************************************************
!!*****
