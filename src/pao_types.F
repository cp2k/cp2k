!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types used by the PAO machinery
!> \author Ole Schuett
! *****************************************************************************
MODULE pao_types
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE linesearch,                      ONLY: linesearch_type
#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_types'

  PUBLIC :: pao_env_type

  TYPE pao_env_type
     ! input values
     REAL(KIND=dp)                    :: check_unitary_tol = 0.0_dp
     REAL(KIND=dp)                    :: check_grad_param_tol = 0.0_dp
     REAL(KIND=dp)                    :: check_grad_full_tol = 0.0_dp
     REAL(KIND=dp)                    :: eps_pao = 0.0_dp
     REAL(KIND=dp)                    :: min_step = 0.0_dp
     REAL(KIND=dp)                    :: mixing = 1.0_dp
     INTEGER                          :: max_pao = -1
     INTEGER                          :: max_cycles = -1
     INTEGER                          :: max_outer_pao = -1
     INTEGER                          :: parameterization = -1
     INTEGER                          :: cg_init_steps = -1
     INTEGER                          :: random_seed = -1
     CHARACTER(LEN=default_string_length) :: preopt_dm_file = ""

     ! output units
     INTEGER                          :: iw = -1
     INTEGER                          :: iw_cg = -1

     ! state variable
     INTEGER                          :: istep = -1
     INTEGER                          :: last_reset = 0
     REAL(KIND=dp)                    :: energy_prev = 0.0_dp
     REAL(KIND=dp)                    :: step_start_time = 0.0_dp
     TYPE(linesearch_type)            :: linesearch

     ! matrices
     TYPE(cp_dbcsr_type)              :: matrix_X
     TYPE(cp_dbcsr_type)              :: matrix_U
     TYPE(cp_dbcsr_type)              :: matrix_U0
     TYPE(cp_dbcsr_type)              :: matrix_Y
     TYPE(cp_dbcsr_type)              :: matrix_N
     TYPE(cp_dbcsr_type)              :: matrix_N_inv
     TYPE(cp_dbcsr_type)              :: matrix_X_orig
     TYPE(cp_dbcsr_type)              :: matrix_G
     TYPE(cp_dbcsr_type)              :: matrix_G_prev
     TYPE(cp_dbcsr_type)              :: matrix_D

     TYPE(cp_dbcsr_type), DIMENSION(:), POINTER :: matrix_V_terms => Null()
  END TYPE

END MODULE pao_types
