!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types used by the PAO machinery
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_types
   USE cp_dbcsr_interface,              ONLY: cp_dbcsr_release,&
                                              cp_dbcsr_type
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE linesearch,                      ONLY: linesearch_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_types'

   PUBLIC :: pao_env_type, training_matrix_type, pao_finalize

   TYPE filename_type
      CHARACTER(LEN=default_path_length) :: fn
   END TYPE filename_type

   TYPE training_matrix_type
      CHARACTER(LEN=default_string_length)      :: kindname = ""
      REAL(dp), DIMENSION(:, :), ALLOCATABLE    :: inputs
      REAL(dp), DIMENSION(:, :), ALLOCATABLE    :: outputs
      REAL(dp), DIMENSION(:), ALLOCATABLE       :: prior
      REAL(dp), DIMENSION(:, :, :), ALLOCATABLE :: NN ! Neural Network
      REAL(dp), DIMENSION(:, :), ALLOCATABLE    :: GP ! Gaussian Process
   END TYPE training_matrix_type

   TYPE pao_env_type
      ! input values
      REAL(KIND=dp)                    :: eps_pao = 0.0_dp
      REAL(KIND=dp)                    :: cg_reset_limit = 0.1_dp
      REAL(KIND=dp)                    :: mixing = 0.0_dp
      REAL(KIND=dp)                    :: regularization = 0.0_dp
      REAL(KIND=dp)                    :: penalty_dist = 0.0_dp
      REAL(KIND=dp)                    :: penalty_strength = 0.0_dp
      REAL(KIND=dp)                    :: check_unitary_tol = 0.0_dp
      REAL(KIND=dp)                    :: check_grad_tol = 0.0_dp
      REAL(KIND=dp)                    :: num_grad_eps = 0.0_dp
      REAL(KIND=dp)                    :: eps_pgf = 0.0_dp
      REAL(KIND=dp)                    :: linpot_precon_delta = 0.0_dp
      REAL(KIND=dp)                    :: linpot_init_delta = 0.0_dp
      REAL(KIND=dp)                    :: linpot_regu_delta = 0.0_dp
      REAL(KIND=dp)                    :: linpot_regu_strength = 0.0_dp
      INTEGER                          :: num_grad_order = -1
      INTEGER                          :: max_pao = -1
      INTEGER                          :: max_cycles = -1
      INTEGER                          :: write_cycles = -1
      INTEGER                          :: parameterization = -1
      INTEGER                          :: optimizer = -1
      INTEGER                          :: cg_init_steps = -1
      CHARACTER(LEN=default_path_length) :: preopt_dm_file = ""
      CHARACTER(LEN=default_path_length) :: restart_file = ""
      TYPE(filename_type), DIMENSION(:), ALLOCATABLE :: ml_training_set

      INTEGER                          :: ml_method = -1
      INTEGER                          :: ml_prior = -1
      INTEGER                          :: ml_descriptor = -1
      REAL(KIND=dp)                    :: ml_tolerance = 0.0_dp
      REAL(KIND=dp)                    :: gp_noise_var = 0.0_dp
      REAL(KIND=dp)                    :: gp_scale = 0.0_dp

      ! output units
      INTEGER                          :: iw = -1
      INTEGER                          :: iw_opt = -1

      ! state variable
      INTEGER                          :: istep = -1
      REAL(KIND=dp)                    :: energy_prev = 0.0_dp
      REAL(KIND=dp)                    :: step_start_time = 0.0_dp
      REAL(KIND=dp)                    :: norm_G = 0.0_dp
      TYPE(linesearch_type)            :: linesearch
      LOGICAL                          :: matrix_X_ready = .FALSE.
      LOGICAL                          :: matrix_P_ready = .FALSE.
      LOGICAL                          :: matrices_YNH0_ready = .FALSE.
      LOGICAL                          :: need_initial_scf = .FALSE.
      LOGICAL                          :: precondition = .FALSE.

      ! matrices
      TYPE(cp_dbcsr_type)              :: matrix_X
      TYPE(cp_dbcsr_type)              :: matrix_U
      TYPE(cp_dbcsr_type)              :: matrix_U0
      TYPE(cp_dbcsr_type)              :: matrix_H0
      TYPE(cp_dbcsr_type)              :: matrix_Y
      TYPE(cp_dbcsr_type)              :: matrix_N
      TYPE(cp_dbcsr_type)              :: matrix_N_inv
      TYPE(cp_dbcsr_type)              :: matrix_X_orig
      TYPE(cp_dbcsr_type)              :: matrix_G
      TYPE(cp_dbcsr_type)              :: matrix_G_prev
      TYPE(cp_dbcsr_type)              :: matrix_D
      TYPE(cp_dbcsr_type)              :: matrix_D_preconed
      TYPE(cp_dbcsr_type)              :: matrix_V_terms
      TYPE(cp_dbcsr_type)              :: matrix_BFGS
      TYPE(cp_dbcsr_type)              :: matrix_precon
      TYPE(cp_dbcsr_type)              :: matrix_precon_inv
      TYPE(cp_dbcsr_type)              :: matrix_R

      TYPE(training_matrix_type), ALLOCATABLE, &
         DIMENSION(:)                  :: ml_training_matrices

   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief Finalize the PAO environment
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_finalize(pao)
      TYPE(pao_env_type)                                 :: pao

      IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| finalize"

      CALL cp_dbcsr_release(pao%matrix_X)
      CALL cp_dbcsr_release(pao%matrix_Y)
      CALL cp_dbcsr_release(pao%matrix_N)
      CALL cp_dbcsr_release(pao%matrix_N_inv)
      CALL cp_dbcsr_release(pao%matrix_H0)

      DEALLOCATE (pao%ml_training_set)
      IF (ALLOCATED(pao%ml_training_matrices)) &
         DEALLOCATE (pao%ml_training_matrices)

      !TODO: should finish printkey

   END SUBROUTINE pao_finalize

END MODULE pao_types
