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
                                              section_vals_type,&
                                              section_vals_val_get
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
                                               pao_gth_param = 104

CONTAINS

! **************************************************************************************************
!> \brief Declare the PAO input section
!> \param pao ...
!> \param input ...
! **************************************************************************************************
   SUBROUTINE parse_pao_section(pao, input)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(section_vals_type), POINTER                   :: input

      CHARACTER(len=*), PARAMETER :: routineN = 'parse_pao_section', &
         routineP = moduleN//':'//routineN

      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)

      ! open main logger
      logger => cp_get_default_logger()
      pao%iw = cp_print_key_unit_nr(logger, input, "DFT%LS_SCF%PAO%PRINT%RUN_INFO", extension=".paolog")
      pao%iw_cg = cp_print_key_unit_nr(logger, input, "DFT%LS_SCF%PAO%PRINT%CG_INFO", extension=".paocglog")

      ! parse input and print
      CALL section_vals_val_get(input, "DFT%LS_SCF%PAO%EPS_PAO", r_val=pao%eps_pao)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "EPS_PAO", pao%eps_pao

      CALL section_vals_val_get(input, "DFT%LS_SCF%PAO%MIXING", r_val=pao%mixing)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "MIXING", pao%mixing

      CALL section_vals_val_get(input, "DFT%LS_SCF%PAO%MAX_PAO", i_val=pao%max_pao)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,I10)") " PAO|", "MAX_PAO", pao%max_pao

      CALL section_vals_val_get(input, "DFT%LS_SCF%PAO%MAX_CYCLES", i_val=pao%max_cycles)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,I10)") " PAO|", "MAX_CYCLES", pao%max_cycles

      CALL section_vals_val_get(input, "DFT%LS_SCF%PAO%PARAMETERIZATION", i_val=pao%parameterization)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,A10)") " PAO|", "PARAMETERIZATION", id2str(pao%parameterization)

      CALL section_vals_val_get(input, "DFT%LS_SCF%PAO%PENALTY_DISTANCE", r_val=pao%penalty_dist)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "PENALTY_DISTANCE", pao%penalty_dist
      IF (pao%penalty_dist < 0.0_dp) CPABORT("PAO: PENALTY_DISTANCE < 0")

      CALL section_vals_val_get(input, "DFT%LS_SCF%PAO%PENALTY_STRENGTH", r_val=pao%penalty_strength)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "PENALTY_STRENGTH", pao%penalty_strength
      IF (pao%penalty_strength < 0.0_dp) CPABORT("PAO: PENALTY_STRENGTH < 0")

      CALL section_vals_val_get(input, "DFT%LS_SCF%PAO%CG_INIT_STEPS", i_val=pao%cg_init_steps)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,I10)") " PAO|", "CG_INIT_STEPS", pao%cg_init_steps
      IF (pao%cg_init_steps < 1) CPABORT("PAO: CG_INIT_STEPS < 1")

      CALL section_vals_val_get(input, "DFT%LS_SCF%PAO%CG_RESET_LIMIT", r_val=pao%cg_reset_limit)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "CG_RESET_LIMIT", pao%cg_reset_limit
      IF (pao%cg_reset_limit < 0.0_dp) CPABORT("PAO: CG_RESET_LIMIT < 0")

      CALL section_vals_val_get(input, "DFT%LS_SCF%PAO%CHECK_UNITARY_TOL", r_val=pao%check_unitary_tol)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "CHECK_UNITARY_TOL", pao%check_unitary_tol

      CALL section_vals_val_get(input, "DFT%LS_SCF%PAO%CHECK_GRADIENT_TOL", r_val=pao%check_grad_tol)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "CHECK_GRADIENT_TOL", pao%check_grad_tol

      CALL section_vals_val_get(input, "DFT%LS_SCF%PAO%NUM_GRADIENT_ORDER", i_val=pao%num_grad_order)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,I10)") " PAO|", "NUM_GRADIENT_ORDER", pao%num_grad_order

      CALL section_vals_val_get(input, "DFT%LS_SCF%PAO%NUM_GRADIENT_EPS", r_val=pao%num_grad_eps)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,T71,E10.1)") " PAO|", "NUM_GRADIENT_EPS", pao%num_grad_eps
      IF (pao%num_grad_eps < 0.0_dp) CPABORT("PAO: NUM_GRADIENT_EPS < 0")

      CALL section_vals_val_get(input, "DFT%LS_SCF%PAO%PREOPT_DM_FILE", c_val=pao%preopt_dm_file)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,A)") " PAO|", "PREOPT_DM_FILE ", pao%preopt_dm_file

      CALL section_vals_val_get(input, "DFT%LS_SCF%PAO%RESTART_FILE", c_val=pao%restart_file)
      IF (pao%iw > 0) WRITE (pao%iw, "(A,T40,A,A)") " PAO|", "RESTART_FILE ", pao%restart_file

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
      CASE DEFAULT
         CPABORT("PAO: unkown id")
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
      TYPE(section_type), POINTER                        :: printkey, subsection

      NULLIFY (keyword, subsection, printkey)

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
                          description="Mixing fraction of new and old solutions after CG", &
                          default_r_val=1.0_dp)
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
                          enum_c_vals=s2a("GTH", "ROTINV", "FOCK", "EXP"), &
                          enum_i_vals=(/pao_gth_param, pao_rotinv_param, pao_fock_param, pao_exp_param/), &
                          enum_desc=s2a("Parametrization based on GTH pseudo potentials", &
                                        "Rotational invariant parametrization", &
                                        "Fock matrix parametrization", &
                                        "Original matrix exponential parametrization"), &
                          default_i_val=pao_rotinv_param)
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

      ! Output *******************************************************************
      CALL section_create(subsection, name="PRINT", &
                          description="Printkey section", &
                          n_keywords=0, n_subsections=1, repeats=.TRUE.)

      CALL cp_print_key_section_create(printkey, "RUN_INFO", &
                                       description="Normal output by PAO", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, printkey)
      CALL section_release(printkey)

      CALL cp_print_key_section_create(printkey, "CG_INFO", &
                                       description="Output by the CG optimizer", &
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

      CALL section_add_subsection(subsection, printkey)
      CALL section_release(printkey)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! CG stuff *****************************************************************
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
