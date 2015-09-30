!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types used by the PAO machinery
!> \author Ole Schuett
! *****************************************************************************
MODULE pao_types
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_type
  USE kinds,                           ONLY: dp
  USE linesearch,                      ONLY: linesearch_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_types'

  PUBLIC :: pao_env_type

  TYPE pao_env_type
     TYPE(cp_dbcsr_type)              :: matrix_X
     TYPE(cp_dbcsr_type)              :: matrix_U
     TYPE(cp_dbcsr_type)              :: matrix_U0
     TYPE(cp_dbcsr_type)              :: matrix_Y
     TYPE(cp_dbcsr_type)              :: matrix_N
     TYPE(cp_dbcsr_type)              :: matrix_N_inv
     TYPE(linesearch_type)            :: linesearch
     TYPE(cp_dbcsr_type)              :: matrix_X_orig
     TYPE(cp_dbcsr_type)              :: matrix_G
     TYPE(cp_dbcsr_type)              :: matrix_G_prev
     TYPE(cp_dbcsr_type)              :: matrix_D
     TYPE(cp_dbcsr_type)              :: matrix_H0
     TYPE(cp_dbcsr_type), DIMENSION(:), POINTER :: matrix_P_guess
     REAL(KIND=dp)                    :: eps_pao
     INTEGER                          :: max_pao
     INTEGER                          :: max_outer_pao
     LOGICAL                          :: line_search_start
     REAL(KIND=dp)                    :: line_search_step
     INTEGER                          :: parameterization
     INTEGER                          :: cg_init_steps
     INTEGER                          :: istep = -1
     INTEGER                          :: iw = -1
     INTEGER                          :: iw_cg = -1
     INTEGER                          :: last_reset = 0
     REAL(KIND=dp)                    :: energy_prev = 0.0_dp
     REAL(KIND=dp)                    :: check_gradient_tol = 0.0_dp
     LOGICAL                          :: check_gradient_abort = .TRUE.
     TYPE(cp_dbcsr_type)              :: matrix_z2c
     TYPE(cp_dbcsr_type), DIMENSION(:), POINTER :: matrix_V_terms => Null()
     TYPE(cp_dbcsr_type)              :: matrix_precon
     TYPE(cp_dbcsr_type)              :: matrix_precon_inv
     INTEGER                          :: preconditioner = -1
  END TYPE

END MODULE pao_types
