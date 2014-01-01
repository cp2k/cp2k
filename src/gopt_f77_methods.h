!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief evaluete the potential energy and its gradients using an array
!>      with same dimension as the particle_set
!> \param gopt_env the geometry optimization environment
!> \param x the position where the function should be evaluated
!> \param f the function value
!> \param gradient the value of its gradient
!> \param error variable to control error logging, stopping,...
!>             see module cp_error_handling
!> \par History
!>      none
!> \author Teodoro Laino [tlaino] - University of Zurich - 01.2008
! *****************************************************************************
INTERFACE

  RECURSIVE SUBROUTINE cp_eval_at(gopt_env,x,f,gradient,master,&
                                  final_evaluation,para_env,error)

    USE cp_para_types,                   ONLY: cp_para_env_type
    USE gopt_f_types,                    ONLY: gopt_f_type
    USE kinds,                           ONLY: dp
#include "cp_common_uses.h"
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x
    REAL(KIND=dp), INTENT(out), OPTIONAL     :: f
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: gradient
    INTEGER, INTENT(IN)                      :: master
    LOGICAL, INTENT(IN), OPTIONAL            :: final_evaluation
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

  END SUBROUTINE cp_eval_at

END INTERFACE
