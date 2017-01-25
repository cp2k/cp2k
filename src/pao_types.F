!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types used by the PAO machinery
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_types
   USE dbcsr_api,                       ONLY: dbcsr_distribution_release,&
                                              dbcsr_distribution_type,&
                                              dbcsr_release,&
                                              dbcsr_type
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

! **************************************************************************************************
!> \brief PAO machine learning data for one atomic kind
!> \var kindname                  name of atomic kind
!> \var inputs                    training points
!> \var outputs                   training points
!> \var prior                     constant prior which is added to prediction
!> \var NN                        trained neural network
!> \var GP                        trained gaussian process
! **************************************************************************************************
   TYPE training_matrix_type
      CHARACTER(LEN=default_string_length)      :: kindname = ""
      REAL(dp), DIMENSION(:, :), ALLOCATABLE    :: inputs
      REAL(dp), DIMENSION(:, :), ALLOCATABLE    :: outputs
      REAL(dp), DIMENSION(:), ALLOCATABLE       :: prior
      REAL(dp), DIMENSION(:, :, :), ALLOCATABLE :: NN ! Neural Network
      REAL(dp), DIMENSION(:, :), ALLOCATABLE    :: GP ! Gaussian Process
   END TYPE training_matrix_type

! **************************************************************************************************
!> \brief The PAO environment type
!> \var eps_pao                   parsed input parameter
!> \var cg_reset_limit            parsed input parameter
!> \var mixing                    parsed input parameter
!> \var regularization            parsed input parameter
!> \var penalty_dist              parsed input parameter
!> \var penalty_strength          parsed input parameter
!> \var check_unitary_tol         parsed input parameter
!> \var check_grad_tol            parsed input parameter
!> \var num_grad_eps              parsed input parameter
!> \var eps_pgf                   parsed input parameter
!> \var linpot_precon_delta       parsed input parameter
!> \var linpot_init_delta         parsed input parameter
!> \var linpot_regu_delta         parsed input parameter
!> \var linpot_regu_strength      parsed input parameter
!> \var num_grad_order            parsed input parameter
!> \var max_pao                   parsed input parameter
!> \var max_cycles                parsed input parameter
!> \var write_cycles              parsed input parameter
!> \var parameterization          parsed input parameter
!> \var optimizer                 parsed input parameter
!> \var cg_init_steps             parsed input parameter
!> \var preopt_dm_file            parsed input parameter
!> \var restart_file              parsed input parameter
!> \var ml_training_set           parsed input parameter
!> \var ml_method                 parsed input parameter
!> \var ml_prior                  parsed input parameter
!> \var ml_descriptor             parsed input parameter
!> \var ml_tolerance              parsed input parameter
!> \var gp_noise_var              parsed input parameter
!> \var gp_scale                  parsed input parameter
!> \var precondition              parsed input parameter
!> \var iw                        output unit for pao in general
!> \var iw_atoms                  output unit for one line summary for each atom
!> \var iw_gap                    output unit for gap of the fock matrix
!> \var iw_fockev                 output unit for eigenvalues of the fock matrix
!> \var iw_opt                    output unit for pao optimizer
!> \var iw_mlvar                  output unit for variances of machine learning predictions
!> \var iw_mldata                 output unit for dumping training data used for machine learning
!> \var istep                     counts pao iterations, ie. number of pao energy evaluations
!> \var energy_prev               energy of previous pao step
!> \var step_start_time           timestamp of when current pao step started
!> \var norm_G                    frobenius-norm of matrix_G or matrix_G_preconed
!> \var linesearch                holds linesearch state
!> \var matrix_X_ready            set when matrix_X is initialized
!> \var matrix_P_ready            set when density matrix is initialized
!> \var constants_ready           set when stuff, which does not depend of atomic positions is ready
!> \var need_initial_scf          set when the initial density matrix is not self-consistend
!> \var matrix_X                  parameters of pao basis, which eventually determine matrix_U. Uses diag_distribution.
!> \var matrix_U                  roation matrix derived from matrix_X. Uses diag_distribution.
!> \var matrix_U0                 constant pre-rotation which serves as initial guess for exp-parametrization. Uses diag_distribution.
!> \var matrix_H0                 Diagonal blocks of core hamiltonian, uses diag_distribution
!> \var matrix_Y                  selector matrix which translates between primary and pao basis.
!>                                basically a block diagonal "rectangular identity matrix". Uses s_matrix-distribution.
!> \var matrix_N                  diagonal matrix filled with 1/sqrt(S) from primary overlap matrix. Uses s_matrix-distribution.
!> \var matrix_N_inv              diagonal matrix filled with sqrt(S) from primary overlap matrix. Uses s_matrix-distribution.
!> \var matrix_N_diag             copy of matrix_N using diag_distribution
!> \var matrix_X_orig             copy made of matrix_X at beginning of optimization cylce, used for mixing. Uses diag_distribution.
!> \var matrix_G                  derivative of pao-energy wrt to matrix_X, ie. the pao-gradient. Uses diag_distribution.
!> \var matrix_G_prev             copy of gradient from previous step, used for conjugate gradient method. Uses diag_distribution.
!> \var matrix_D                  Current line-search direction,  used for conjugate gradient method. Uses diag_distribution.
!> \var matrix_D_preconed         Current line-search direction with preconditioner applied.
!>                                This copy is keept, because application of inverse preconditioner
!>                                introduces too much numeric noise. Uses diag_distribution.
!> \var matrix_V_terms            Potential terms, used by linpot and gth parametrization, Uses diag_distribution.
!> \var matrix_BFGS               Approximate inverse hessian, used by BFGS optimizer, Uses diag_distribution.
!> \var matrix_precon             preconditioner, uses diag_distribution.
!> \var matrix_precon_inv         inverse of matrix_precon, uses diag_distribution.
!> \var matrix_R                  Rgularization, uses diag_distribution
!> \var ml_training_matrices      holds training data and trained machine learning model
!> \var diag_distribution         DBCSR distribution to spreads diagonal blocks evenly across ranks
! **************************************************************************************************
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
      LOGICAL                          :: precondition = .FALSE.
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
      INTEGER                          :: iw_atoms = -1
      INTEGER                          :: iw_gap = -1
      INTEGER                          :: iw_fockev = -1
      INTEGER                          :: iw_opt = -1
      INTEGER                          :: iw_mlvar = -1
      INTEGER                          :: iw_mldata = -1

      ! state variable
      INTEGER                          :: istep = -1
      REAL(KIND=dp)                    :: energy_prev = 0.0_dp
      REAL(KIND=dp)                    :: step_start_time = 0.0_dp
      REAL(KIND=dp)                    :: norm_G = 0.0_dp
      TYPE(linesearch_type)            :: linesearch
      LOGICAL                          :: matrix_X_ready = .FALSE.
      LOGICAL                          :: matrix_P_ready = .FALSE.
      LOGICAL                          :: constants_ready = .FALSE.
      LOGICAL                          :: need_initial_scf = .FALSE.

      ! matrices
      TYPE(dbcsr_type)              :: matrix_X
      TYPE(dbcsr_type)              :: matrix_U
      TYPE(dbcsr_type)              :: matrix_U0
      TYPE(dbcsr_type)              :: matrix_H0
      TYPE(dbcsr_type)              :: matrix_Y
      TYPE(dbcsr_type)              :: matrix_N
      TYPE(dbcsr_type)              :: matrix_N_inv
      TYPE(dbcsr_type)              :: matrix_N_diag
      TYPE(dbcsr_type)              :: matrix_X_orig
      TYPE(dbcsr_type)              :: matrix_G
      TYPE(dbcsr_type)              :: matrix_G_prev
      TYPE(dbcsr_type)              :: matrix_D
      TYPE(dbcsr_type)              :: matrix_D_preconed
      TYPE(dbcsr_type)              :: matrix_V_terms
      TYPE(dbcsr_type)              :: matrix_BFGS
      TYPE(dbcsr_type)              :: matrix_precon
      TYPE(dbcsr_type)              :: matrix_precon_inv
      TYPE(dbcsr_type)              :: matrix_R

      TYPE(training_matrix_type), ALLOCATABLE, &
         DIMENSION(:)                  :: ml_training_matrices

      TYPE(dbcsr_distribution_type)     :: diag_distribution
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief Finalize the PAO environment
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_finalize(pao)
      TYPE(pao_env_type)                                 :: pao

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_finalize', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL dbcsr_release(pao%matrix_X)
      CALL dbcsr_release(pao%matrix_Y)
      CALL dbcsr_release(pao%matrix_N)
      CALL dbcsr_release(pao%matrix_N_inv)
      CALL dbcsr_release(pao%matrix_N_diag)
      CALL dbcsr_release(pao%matrix_H0)

      DEALLOCATE (pao%ml_training_set)
      IF (ALLOCATED(pao%ml_training_matrices)) &
         DEALLOCATE (pao%ml_training_matrices)

      CALL dbcsr_distribution_release(pao%diag_distribution)

      !TODO: should finish printkey

      CALL timestop(handle)

   END SUBROUTINE pao_finalize

END MODULE pao_types
