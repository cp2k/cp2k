!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE fit_parameters_types

  USE kinds,                           ONLY: dp

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: linklist_real_type, linklist_parameters_type, &
            fit_parameter_type, fit_info_type, fit_averages_type

  TYPE linklist_real_type
     REAL (KIND=dp), POINTER :: real_value
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
     REAL (KIND=dp) :: max_tol
  END TYPE fit_parameter_type
 
  TYPE fit_info_type
     REAL (KIND=dp), DIMENSION ( 3 ) :: f_ref, f_fit, r
  END TYPE fit_info_type

  TYPE fit_averages_type
    REAL(KIND=dp), POINTER, DIMENSION(:,:) :: ave_sig, ave_eps
    REAL(KIND=dp), POINTER, DIMENSION(:,:) :: ave_a, ave_b, ave_c
  END TYPE fit_averages_type
                                                                                          

END MODULE fit_parameters_types

!******************************************************************************
!!*****
