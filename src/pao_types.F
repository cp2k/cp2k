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
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE linesearch,                      ONLY: linesearch_type
#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_types'

  PUBLIC :: pao_env_type
  PUBLIC :: mcweeny_history_type, mcweeny_history_p_type

  TYPE mcweeny_history_type
    TYPE(cp_dbcsr_type)                        :: m
    INTEGER                                    :: count = -1
    TYPE(mcweeny_history_type), POINTER        :: next => Null()
  END TYPE mcweeny_history_type

  TYPE mcweeny_history_p_type
    TYPE(mcweeny_history_type), POINTER        :: p => Null()
  END TYPE mcweeny_history_p_type

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
     TYPE(cp_dbcsr_type), DIMENSION(:), POINTER :: matrix_P_guess => Null()
     REAL(KIND=dp)                    :: eps_pao = 0.0_dp
     REAL(KIND=dp)                    :: eps_preopt = 0.0_dp
     INTEGER                          :: max_pao = -1
     INTEGER                          :: max_outer_pao = -1
     INTEGER                          :: parameterization = -1
     INTEGER                          :: cg_init_steps = -1
     INTEGER                          :: istep = -1
     INTEGER                          :: iw = -1
     INTEGER                          :: iw_cg = -1
     INTEGER                          :: last_reset = 0
     REAL(KIND=dp)                    :: energy_prev = 0.0_dp
     REAL(KIND=dp)                    :: check_unitary_tol = 0.0_dp
     REAL(KIND=dp)                    :: check_grad_param_tol = 0.0_dp
     REAL(KIND=dp)                    :: check_grad_mcweeny_tol = 0.0_dp
     TYPE(cp_dbcsr_type)              :: matrix_z2c
     TYPE(cp_dbcsr_type), DIMENSION(:), POINTER :: matrix_V_terms => Null()
     TYPE(cp_dbcsr_type)              :: matrix_precon
     TYPE(cp_dbcsr_type)              :: matrix_precon_inv
     INTEGER                          :: preconditioner = -1
     CHARACTER(LEN=default_string_length) :: preopt_dm_file = ""
     TYPE(mcweeny_history_p_type),&
      DIMENSION(:), POINTER           :: mcweeny_history => Null()
     REAL(KIND=dp)                    :: step_start_time = 0.0_dp
  END TYPE

END MODULE pao_types
