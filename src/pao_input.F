!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE pao_input
   USE bibliography,                    ONLY: Berghold2011
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: add_last_numeric,&
                                              cp_print_key_section_create,&
                                              cp_print_key_unit_nr,&
                                              high_print_level,&
                                              low_print_level
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_add_subsection,&
                                              section_create,&
                                              section_release,&
                                              section_type,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE input_val_types,                 ONLY: lchar_t,&
                                              real_t
   USE kinds,                           ONLY: dp
   USE linesearch,                      ONLY: linesearch_create_section
   USE pao_types,                       ONLY: pao_env_type
   USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_input'

   PUBLIC :: create_pao_section, parse_pao_section, id2str

   INTEGER, PARAMETER, PUBLIC               :: pao_rotinv_param = 101, &
                                               pao_fock_param = 102, &
                                               pao_exp_param = 103, &
                                               pao_gth_param = 104, &
                                               pao_opt_cg = 301, &
                                               pao_opt_bfgs = 302, &
                                               pao_ml_gp = 401, &
                                               pao_ml_nn = 402, &
                                               pao_ml_lazy = 403, &
                                               pao_ml_prior_zero = 501, &
                                               pao_ml_prior_mean = 502, &
                                               pao_ml_desc_pot = 601, &
                                               pao_ml_desc_overlap = 602, &
                                               pao_ml_desc_r12 = 603

CONTAINS

! **************************************************************************************************
!> \brief Declare the PAO input section
!> \param pao ...
!> \param input ...
! **************************************************************************************************
   SUBROUTINE parse_pao_section(pao, input)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(section_vals_type), POINTER                   :: input

      INTEGER                                            :: i, n_rep, ntrainfiles
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: pao_section, training_set_section

      pao_section => section_vals_get_subs_vals(input, "DFT%LS_SCF%PAO")

      ! open main logger
      logger => cp_get_default_logger()
      pao%iw = cp_print_key_unit_nr(logger, pao_section, "PRINT%RUN_INFO", extension=".paolog")
      pao%iw_atoms = cp_print_key_unit_nr(logger, pao_section, "PRINT%ATOM_INFO", extension=".paolog")
      pao%iw_gap = cp_print_key_unit_nr(logger, pao_section, "PRINT%FOCK_GAP", extension=".paolog")
      pao%iw_fockev = cp_print_key_unit_nr(logger, pao_section, "PRINT%FOCK_EIGENVALUES", extension=".paolog")
      pao%iw_opt = cp_print_key_unit_nr(logger, pao_section, "PRINT%OPT_INFO", extension=".paolog")
      pao%iw_mlvar = cp_print_key_unit_nr(logger, pao_section, "PRINT%ML_VARIANCE", extension=".paolog")
      pao%iw_mldata = cp_print_key_unit_nr(logger, pao_section, "PRINT%ML_TRAINING_DATA", extension=".paolog")

      IF (pao%iw > 0) WRITE (pao%iw, *) "" ! an empty separator line

      ! parse input and print

      CALL section_vals_val_get(pao_section, "EPS_PAO", r_val=pao%eps_pao)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "EPS_PAO", pao%eps_pao

      CALL section_vals_val_get(pao_section, "MIXING", r_val=pao%mixing)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "MIXING", pao%mixing

      CALL section_vals_val_get(pao_section, "MAX_PAO", i_val=pao%max_pao)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,I10)") " PAO|", "MAX_PAO", pao%max_pao

      CALL section_vals_val_get(pao_section, "MAX_CYCLES", i_val=pao%max_cycles)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,I10)") " PAO|", "MAX_CYCLES", pao%max_cycles

      CALL section_vals_val_get(pao_section, "PARAMETERIZATION", i_val=pao%parameterization)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,A10)") " PAO|", "PARAMETERIZATION", id2str(pao%parameterization)

      CALL section_vals_val_get(pao_section, "PRECONDITION", l_val=pao%precondition)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,L10)") " PAO|", "PRECONDITION", pao%precondition

      CALL section_vals_val_get(pao_section, "REGULARIZATION", r_val=pao%regularization)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "REGULARIZATION", pao%regularization
      IF (pao%regularization < 0.0_dp) CPABORT("PAO: REGULARIZATION < 0")

      CALL section_vals_val_get(input, "DFT%QS%EPS_DEFAULT", r_val=pao%eps_pgf) ! default value
      CALL section_vals_val_get(pao_section, "EPS_PGF", n_rep_val=n_rep)
      IF (n_rep /= 0) CALL section_vals_val_get(pao_section, "EPS_PGF", r_val=pao%eps_pgf)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "EPS_PGF", pao%eps_pgf
      IF (pao%eps_pgf < 0.0_dp) CPABORT("PAO: EPS_PGF < 0")

      CALL section_vals_val_get(pao_section, "PENALTY_DISTANCE", r_val=pao%penalty_dist)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "PENALTY_DISTANCE", pao%penalty_dist
      IF (pao%penalty_dist < 0.0_dp) CPABORT("PAO: PENALTY_DISTANCE < 0")

      CALL section_vals_val_get(pao_section, "PENALTY_STRENGTH", r_val=pao%penalty_strength)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "PENALTY_STRENGTH", pao%penalty_strength
      IF (pao%penalty_strength < 0.0_dp) CPABORT("PAO: PENALTY_STRENGTH < 0")

      CALL section_vals_val_get(pao_section, "LINPOT_PRECONDITION_DELTA", r_val=pao%linpot_precon_delta)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "LINPOT_PRECONDITION_DELTA", pao%linpot_precon_delta
      IF (pao%linpot_precon_delta < 0.0_dp) CPABORT("PAO: LINPOT_PRECONDITION_DELTA < 0")

      CALL section_vals_val_get(pao_section, "LINPOT_INITGUESS_DELTA", r_val=pao%linpot_init_delta)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "LINPOT_INITGUESS_DELT", pao%linpot_init_delta
      IF (pao%linpot_init_delta < 0.0_dp) CPABORT("PAO: LINPOT_INITGUESS_DELTA < 0")

      CALL section_vals_val_get(pao_section, "LINPOT_REGULARIZATION_DELTA", r_val=pao%linpot_regu_delta)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "LINPOT_REGULARIZATION_DELTA", pao%linpot_regu_delta
      IF (pao%linpot_regu_delta < 0.0_dp) CPABORT("PAO: LINPOT_REGULARIZATION_DELTA < 0")

      CALL section_vals_val_get(pao_section, "LINPOT_REGULARIZATION_STRENGTH", r_val=pao%linpot_regu_strength)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "LINPOT_REGULARIZATION_STRENGTH", pao%linpot_regu_strength
      IF (pao%linpot_regu_strength < 0.0_dp) CPABORT("PAO: LINPOT_REGULARIZATION_STRENGTH < 0")

      CALL section_vals_val_get(pao_section, "OPTIMIZER", i_val=pao%optimizer)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,A10)") " PAO|", "OPTIMIZER", id2str(pao%optimizer)

      CALL section_vals_val_get(pao_section, "CG_INIT_STEPS", i_val=pao%cg_init_steps)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,I10)") " PAO|", "CG_INIT_STEPS", pao%cg_init_steps
      IF (pao%cg_init_steps < 1) CPABORT("PAO: CG_INIT_STEPS < 1")

      CALL section_vals_val_get(pao_section, "CG_RESET_LIMIT", r_val=pao%cg_reset_limit)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "CG_RESET_LIMIT", pao%cg_reset_limit
      IF (pao%cg_reset_limit < 0.0_dp) CPABORT("PAO: CG_RESET_LIMIT < 0")

      CALL section_vals_val_get(pao_section, "CHECK_UNITARY_TOL", r_val=pao%check_unitary_tol)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "CHECK_UNITARY_TOL", pao%check_unitary_tol

      CALL section_vals_val_get(pao_section, "CHECK_GRADIENT_TOL", r_val=pao%check_grad_tol)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "CHECK_GRADIENT_TOL", pao%check_grad_tol

      CALL section_vals_val_get(pao_section, "NUM_GRADIENT_ORDER", i_val=pao%num_grad_order)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,I10)") " PAO|", "NUM_GRADIENT_ORDER", pao%num_grad_order

      CALL section_vals_val_get(pao_section, "NUM_GRADIENT_EPS", r_val=pao%num_grad_eps)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "NUM_GRADIENT_EPS", pao%num_grad_eps
      IF (pao%num_grad_eps < 0.0_dp) CPABORT("PAO: NUM_GRADIENT_EPS < 0")

      CALL section_vals_val_get(pao_section, "PRINT%RESTART%WRITE_CYCLES", i_val=pao%write_cycles)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,I10)") " PAO|", "PRINT%RESTART%WRITE_CYCLES", pao%write_cycles

      CALL section_vals_val_get(pao_section, "RESTART_FILE", c_val=pao%restart_file)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,A)") " PAO|", "RESTART_FILE ", TRIM(pao%restart_file)

      CALL section_vals_val_get(pao_section, "PREOPT_DM_FILE", c_val=pao%preopt_dm_file)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,A)") " PAO|", "PREOPT_DM_FILE ", TRIM(pao%preopt_dm_file)

      CALL section_vals_val_get(pao_section, "MACHINE_LEARNING%METHOD", i_val=pao%ml_method)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,A10)") " PAO|", "MACHINE_LEARNING%METHOD", id2str(pao%ml_method)

      CALL section_vals_val_get(pao_section, "MACHINE_LEARNING%PRIOR", i_val=pao%ml_prior)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,A10)") " PAO|", "MACHINE_LEARNING%PRIOR", id2str(pao%ml_prior)

      CALL section_vals_val_get(pao_section, "MACHINE_LEARNING%DESCRIPTOR", i_val=pao%ml_descriptor)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,A10)") " PAO|", "MACHINE_LEARNING%DESCRIPTOR", id2str(pao%ml_descriptor)

      CALL section_vals_val_get(pao_section, "MACHINE_LEARNING%TOLERANCE", r_val=pao%ml_tolerance)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "MACHINE_LEARNING%TOLERANCE", pao%ml_tolerance

      CALL section_vals_val_get(pao_section, "MACHINE_LEARNING%GP_NOISE_VAR", r_val=pao%gp_noise_var)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "MACHINE_LEARNING%GP_NOISE_VAR", pao%gp_noise_var

      CALL section_vals_val_get(pao_section, "MACHINE_LEARNING%GP_SCALE", r_val=pao%gp_scale)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "MACHINE_LEARNING%GP_SCALE", pao%gp_scale

      ! parse MACHINE_LEARNING%TRAINING_SET section
      training_set_section => section_vals_get_subs_vals(pao_section, "MACHINE_LEARNING%TRAINING_SET")
      CALL section_vals_val_get(training_set_section, "_DEFAULT_KEYWORD_", n_rep_val=ntrainfiles)
      ALLOCATE (pao%ml_training_set(ntrainfiles))
      DO i = 1, ntrainfiles
         CALL section_vals_val_get(training_set_section, "_DEFAULT_KEYWORD_", &
                                   i_rep_val=i, c_val=pao%ml_training_set(i)%fn)
         IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,A)") " PAO|", "MACHINE_LEARNING%TRAINING_SET", &
            TRIM(pao%ml_training_set(i)%fn)
      ENDDO

      IF (pao%iw > 0) WRITE (pao%iw, *) "" ! an empty separator line

   END SUBROUTINE parse_pao_section

! **************************************************************************************************
!> \brief Helper routine
!> \param id ...
!> \retval s ...
! **************************************************************************************************
   FUNCTION id2str(id) RESULT(s)
      INTEGER                                            :: id
      CHARACTER(LEN=10)                                  :: s

      CHARACTER(len=*), PARAMETER :: routineN = 'id2str', routineP = moduleN//':'//routineN

      SELECT CASE (id)
      CASE (pao_gth_param)
         s = "GTH"
      CASE (pao_rotinv_param)
         s = "ROTINV"
      CASE (pao_fock_param)
         s = "FOCK"
      CASE (pao_exp_param)
         s = "EXP"
      CASE (pao_opt_cg)
         s = "CG"
      CASE (pao_opt_bfgs)
         s = "BFGS"
      CASE (pao_ml_gp)
         s = "GAUSSIAN_PROCESS"
      CASE (pao_ml_nn)
         s = "NEURAL_NETWORK"
      CASE (pao_ml_lazy)
         s = "LAZY"
      CASE (pao_ml_prior_zero)
         s = "ZERO"
      CASE (pao_ml_prior_mean)
         s = "MEAN"
      CASE (pao_ml_desc_pot)
         s = "POTENTIAL"
      CASE (pao_ml_desc_overlap)
         s = "OVERLAP"
      CASE (pao_ml_desc_r12)
         s = "R12"
      CASE DEFAULT
         CPABORT("PAO: unknown id")
      END SELECT
      s = ADJUSTR(s)
   END FUNCTION id2str

! **************************************************************************************************
!> \brief Creates the PAO subsection of the linear scaling section.
!> \param section ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE create_pao_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_pao_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: printkey, subsection, subsubsection

      NULLIFY (keyword, subsection, subsubsection, printkey)

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, "PAO", repeats=.FALSE., &
                          description="Polarized Atomic Orbital Method", &
                          citations=(/Berghold2011/))

      ! Convergence Criteria *****************************************************
      CALL keyword_create(keyword, name="EPS_PAO", &
                          description="Convergence criteria for PAO optimization.", &
                          default_r_val=1.e-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="MIXING", &
                          description="Mixing fraction of new and old optimizied solutions.", &
                          default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="MAX_PAO", &
                          description="Maximum number of PAO basis optimization steps.", &
                          default_i_val=1000)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="MAX_CYCLES", &
                          description="Maximum number of PAO line search cycles for a given hamiltonian.", &
                          default_i_val=1000)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Parametrization **********************************************************
      CALL keyword_create(keyword, name="PARAMETERIZATION", &
                          description="Parametrization of the mapping between the primary and the PAO basis.", &
                          enum_c_vals=s2a("ROTINV", "FOCK", "GTH", "EXP"), &
                          enum_i_vals=(/pao_rotinv_param, pao_fock_param, pao_gth_param, pao_exp_param/), &
                          enum_desc=s2a("Rotational invariant parametrization (machine learnable)", &
                                        "Fock matrix parametrization", &
                                        "Parametrization based on GTH pseudo potentials", &
                                        "Original matrix exponential parametrization"), &
                          default_i_val=pao_rotinv_param)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="REGULARIZATION", &
                          description="Strength of regularization term which ensures parameters remain small.", &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="PENALTY_DISTANCE", &
                          description="Distance at which approaching eigenvalues are penalized to prevent degeneration.", &
                          default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="PENALTY_STRENGTH", &
                          description="Strength of the penalty term which prevents degenerate eigenvalues.", &
                          default_r_val=0.005_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="PRECONDITION", &
                          description="Apply a preconditioner to the parametrization.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="EPS_PGF", &
                          description="Sets precision for potential and descriptor matrix elements. "// &
                          "Overrides DFT/QS/EPS_DEFAULT value.", type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Preopt  ******************************************************************
      CALL keyword_create(keyword, name="PREOPT_DM_FILE", &
                          description="Read pre-optimized density matrix from given file.", &
                          repeats=.FALSE., default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Misc ********************************************************************
      CALL keyword_create(keyword, name="RESTART_FILE", &
                          description="Reads given files as restart for PAO basis", &
                          repeats=.FALSE., default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="CHECK_GRADIENT_TOL", &
                          description="Tolerance for check of full analytic gradient against the numeric one."// &
                          " Negative values mean don't check at all.", &
                          default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="NUM_GRADIENT_EPS", &
                          description="Step length used for the numeric derivative when checking the gradient.", &
                          default_r_val=1e-8_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="NUM_GRADIENT_ORDER", &
                          description="Order of the numeric derivative when checking the gradient. "// &
                          "Possible values are 2, 4, and 6.", &
                          default_i_val=2)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="CHECK_UNITARY_TOL", &
                          description="Check if rotation matrix is unitary."// &
                          " Negative values mean don't check at all.", &
                          default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Linpot settings  *********************************************************
      CALL keyword_create(keyword, name="LINPOT_PRECONDITION_DELTA", &
                          description="Eigenvalue threshold used for preconditioning.", &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="LINPOT_INITGUESS_DELTA", &
                          description="Eigenvalue threshold used for calculating initial guess.", &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="LINPOT_REGULARIZATION_DELTA", &
                          description="Eigenvalue threshold used for regularization.", &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="LINPOT_REGULARIZATION_STRENGTH", &
                          description="Strength of regularization on linpot layer.", &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Machine Learning *********************************************************
      CALL section_create(subsection, name="MACHINE_LEARNING", description="Machine learning section")

      CALL keyword_create(keyword, name="METHOD", &
                          description="Machine learning scheme used to predict PAO basis sets.", &
                          enum_c_vals=s2a("GAUSSIAN_PROCESS", "NEURAL_NETWORK", "LAZY"), &
                          enum_i_vals=(/pao_ml_gp, pao_ml_nn, pao_ml_lazy/), &
                          enum_desc=s2a("Gaussian Process", "Neural Network", "Just rely on prior"), &
                          default_i_val=pao_ml_gp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="PRIOR", &
                          description="Prior used for predictions.", &
                          enum_c_vals=s2a("ZERO", "MEAN"), &
                          enum_i_vals=(/pao_ml_prior_zero, pao_ml_prior_mean/), &
                          enum_desc=s2a("Simply use zero", "Use average of training-set"), &
                          default_i_val=pao_ml_prior_zero)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="DESCRIPTOR", &
                          description="Descriptor used as input for machine learning.", &
                          enum_c_vals=s2a("POTENTIAL", "OVERLAP", "R12"), &
                          enum_i_vals=(/pao_ml_desc_pot, pao_ml_desc_overlap, pao_ml_desc_r12/), &
                          enum_desc=s2a("Eigenvalues of local potential matrix", &
                                        "Eigenvalues of local overlap matrix", &
                                        "Distance between two atoms (just for testing)"), &
                          default_i_val=pao_ml_desc_pot)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="TOLERANCE", &
                          description="Maximum variance tolerated when making predictions.", &
                          default_r_val=1.0E-2_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="GP_NOISE_VAR", &
                          description="Variance of noise used for Gaussian Process machine learning.", &
                          default_r_val=0.1_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="GP_SCALE", &
                          description="Length scale used for Gaussian Process machine learning.", &
                          default_r_val=0.05_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      ! special free-text section similar to SUBSYS%COORD
      CALL section_create(subsubsection, name="TRAINING_SET", &
                          description="Lists PAO-restart file used for training")
      CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_", &
                          description="One file name per line.", &
                          repeats=.TRUE., type_of_var=lchar_t)
      CALL section_add_keyword(subsubsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(subsection, subsubsection)
      CALL section_release(subsubsection)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! Output *******************************************************************
      CALL section_create(subsection, name="PRINT", &
                          description="Printkey section", &
                          n_keywords=0, n_subsections=1, repeats=.TRUE.)

      CALL cp_print_key_section_create(printkey, "RUN_INFO", &
                                       description="Normal output by PAO", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, printkey)
      CALL section_release(printkey)

      CALL cp_print_key_section_create(printkey, "ATOM_INFO", &
                                       description="One line summary for each atom", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, printkey)
      CALL section_release(printkey)

      CALL cp_print_key_section_create(printkey, "FOCK_GAP", &
                                       description="Gap of the fock matrix for each atom", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, printkey)
      CALL section_release(printkey)

      CALL cp_print_key_section_create(printkey, "FOCK_EIGENVALUES", &
                                       description="Eigenvalues of the fock matrix for each atom", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, printkey)
      CALL section_release(printkey)

      CALL cp_print_key_section_create(printkey, "ML_VARIANCE", &
                                       description="Variances of machine learning predictions for each atom", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, printkey)
      CALL section_release(printkey)

      CALL cp_print_key_section_create(printkey, "ML_TRAINING_DATA", &
                                       description="Dumps training data used for machine learning", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, printkey)
      CALL section_release(printkey)

      CALL cp_print_key_section_create(printkey, "OPT_INFO", &
                                       description="Output by the optimizer", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, printkey)
      CALL section_release(printkey)

      CALL cp_print_key_section_create(printkey, "RESTART", &
                                       description="Restart file of PAO basis", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="")

      CALL keyword_create(keyword, name="BACKUP_COPIES", &
                          description="Specifies the maximum number of backup copies.", &
                          usage="BACKUP_COPIES {int}", &
                          default_i_val=1)
      CALL section_add_keyword(printkey, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="WRITE_CYCLES", &
                          description="Maximum number of PAO line search cycles until a restart is written.", &
                          default_i_val=100)
      CALL section_add_keyword(printkey, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(subsection, printkey)
      CALL section_release(printkey)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! OPT stuff ****************************************************************
      CALL keyword_create(keyword, name="OPTIMIZER", &
                          description="Optimizer used to find PAO basis.", &
                          enum_c_vals=s2a("CG", "BFGS"), &
                          enum_i_vals=(/pao_opt_cg, pao_opt_bfgs/), &
                          enum_desc=s2a("Conjugate gradient algorithm", &
                                        "Broyden-Fletcher-Goldfarb-Shanno algorithm"), &
                          default_i_val=pao_opt_cg)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="CG_INIT_STEPS", &
                          description="Number of steepest descent steps before starting the"// &
                          " conjugate gradients optimization.", &
                          default_i_val=2)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="CG_RESET_LIMIT", &
                          description="The CG is reseted if the cosine of the angle between the last "// &
                          "search direction and the new gradient is larger that the limit.", &
                          default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL linesearch_create_section(subsection)
      CALL section_add_subsection(section, subsection)

      CALL section_release(subsection)
   END SUBROUTINE create_pao_section

END MODULE pao_input
