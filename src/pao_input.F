!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

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

  INTEGER, PARAMETER, PUBLIC               :: pao_none_precon=421,&
                                              pao_sinv_precon=422

  INTEGER, PARAMETER, PUBLIC               :: pao_rotinv_param=101,&
                                              pao_fock_param=102,&
                                              pao_exp_param=103

CONTAINS

! *****************************************************************************
!> \brief Declare the PAO input section
!> \param pao ...
!> \param input ...
! *****************************************************************************
  SUBROUTINE parse_pao_section(pao, input)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(section_vals_type), POINTER         :: input

    CHARACTER(len=*), PARAMETER :: routineN = 'parse_pao_section', &
      routineP = moduleN//':'//routineN

    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)

    ! open main logger
    logger => cp_get_default_logger()
    pao%iw = cp_print_key_unit_nr(logger, input, "DFT%LS_SCF%PAO%PRINT%RUN_INFO", extension=".paolog")
    pao%iw_cg = cp_print_key_unit_nr(logger, input, "DFT%LS_SCF%PAO%PRINT%CG_INFO", extension=".paocglog")

    ! parse input and print
    CALL section_vals_val_get(input,"DFT%LS_SCF%PAO%EPS_PAO",r_val=pao%eps_pao)
    IF(pao%iw>0) WRITE (pao%iw,"(A,T40,A,T71,E10.1)") " PAO|","EPS_PAO", pao%eps_pao

    CALL section_vals_val_get(input,"DFT%LS_SCF%PAO%MAX_PAO",i_val=pao%max_pao)
    IF(pao%iw>0) WRITE (pao%iw,"(A,T40,A,T71,I10)") " PAO|","MAX_PAO", pao%max_pao

    CALL section_vals_val_get(input,"DFT%LS_SCF%PAO%MAX_OUTER_PAO",i_val=pao%max_outer_pao)
    IF(pao%iw>0) WRITE (pao%iw,"(A,T40,A,T71,I10)") " PAO|","MAX_OUTER_PAO", pao%max_outer_pao

    CALL section_vals_val_get(input,"DFT%LS_SCF%PAO%PARAMETERIZATION",i_val=pao%parameterization)
    IF(pao%iw>0) WRITE (pao%iw,"(A,T40,A,T71,A10)") " PAO|","PARAMETERIZATION", id2str(pao%parameterization)

    CALL section_vals_val_get(input,"DFT%LS_SCF%PAO%PRECONDITIONER",i_val=pao%preconditioner)
    IF(pao%iw>0) WRITE (pao%iw,"(A,T40,A,T71,A10)") " PAO|","PRECONDITIONER", id2str(pao%preconditioner)

    CALL section_vals_val_get(input,"DFT%LS_SCF%PAO%CG_INIT_STEPS",i_val=pao%cg_init_steps)
    IF(pao%iw>0) WRITE (pao%iw,"(A,T40,A,T71,I10)") " PAO|","CG_INIT_STEPS",pao%cg_init_steps
    IF(pao%cg_init_steps<1) CPABORT("PAO: CG_INIT_STEPS < 1")

    CALL section_vals_val_get(input,"DFT%LS_SCF%PAO%CHECK_UNITARY_TOL",r_val=pao%check_unitary_tol)
    IF(pao%iw>0) WRITE (pao%iw,"(A,T40,A,T71,E10.1)") " PAO|","CHECK_UNITARY_TOL",pao%check_unitary_tol

    CALL section_vals_val_get(input,"DFT%LS_SCF%PAO%CHECK_GRADIENT_PARAM_TOL",r_val=pao%check_grad_param_tol)
    IF(pao%iw>0) WRITE (pao%iw,"(A,T40,A,T71,E10.1)") " PAO|","CHECK_GRADIENT_PARAM_TOL",pao%check_grad_param_tol

    CALL section_vals_val_get(input,"DFT%LS_SCF%PAO%CHECK_GRADIENT_MCWEENY_TOL",r_val=pao%check_grad_mcweeny_tol)
    IF(pao%iw>0) WRITE (pao%iw,"(A,T40,A,T71,E10.1)") " PAO|","CHECK_GRADIENT_MCWEENY_TOL",pao%check_grad_mcweeny_tol

    CALL section_vals_val_get(input,"DFT%LS_SCF%PAO%PREOPT_DM_FILE",c_val=pao%preopt_dm_file)
    IF(pao%iw>0) WRITE (pao%iw,"(A,T40,A,A)") " PAO|","PREOPT_DM_FILE ",pao%preopt_dm_file

    CALL section_vals_val_get(input,"DFT%LS_SCF%PAO%EPS_PREOPT",r_val=pao%eps_preopt)
    IF(pao%iw>0) WRITE (pao%iw,"(A,T40,A,T71,E10.1)") " PAO|","EPS_PREOPT", pao%eps_preopt

  END SUBROUTINE parse_pao_section

! *****************************************************************************
!> \brief Helper routine
!> \param id ...
!> \retval s ...
! *****************************************************************************
  FUNCTION id2str(id) RESULT(s)
    INTEGER                                  :: id
    CHARACTER(LEN=10)                        :: s

    CHARACTER(len=*), PARAMETER :: routineN = 'id2str', &
      routineP = moduleN//':'//routineN

    SELECT CASE(id)
      CASE(pao_none_precon)
        s = "NONE"
      CASE(pao_sinv_precon)
        s = "SINV"
      CASE(pao_rotinv_param)
        s = "ROTINV"
      CASE(pao_fock_param)
        s = "FOCK"
      CASE(pao_exp_param)
        s = "EXP"
      CASE DEFAULT
        CPABORT("PAO: unkown id")
    END SELECT
    s = ADJUSTR(s)
  END FUNCTION id2str

! *****************************************************************************
!> \brief Creates the PAO subsection of the linear scaling section.
!> \param section ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE create_pao_section(section)
    TYPE(section_type), POINTER              :: section

    CHARACTER(len=*), PARAMETER :: routineN = 'create_pao_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: printkey, subsection

    failure=.FALSE.
    NULLIFY (keyword, subsection, printkey)

    CPASSERT(.NOT.ASSOCIATED(section))
    CALL section_create(section,"PAO",repeats=.FALSE.,&
                        description="Polarize Atomic Orbital Method",&
                        citations=(/Berghold2011/))

    ! Convergence Criteria *****************************************************
    CALL keyword_create(keyword, name="EPS_PAO",&
         description="Convergence criteria for PAO optimization.",&
         default_r_val=1.e-5_dp)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="MAX_PAO",&
         description="Maximum number of PAO basis optimization steps.",&
         default_i_val=1000)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="MAX_OUTER_PAO",&
         description="Number of steps after which the CG is reseted in any case.",&
         default_i_val=150)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    ! Parametrization **********************************************************
    CALL keyword_create(keyword, name="PRECONDITIONER",&
         description="Preconditioner used by the linear potential parametrization.",&
         enum_c_vals=s2a("NONE","SINV"),&
         enum_i_vals=(/pao_none_precon, pao_sinv_precon/),&
         enum_desc=s2a("No preconditioner",&
                       "Invert overlap matrix between V terms"),&
         default_i_val=pao_sinv_precon)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="PARAMETERIZATION",&
         description="Parametrization of the mapping between the primary and the PAO basis.",&
         enum_c_vals=s2a("ROTINV","FOCK","EXP"),&
         enum_i_vals=(/pao_rotinv_param, pao_fock_param, pao_exp_param/),&
         enum_desc=s2a("Rotational invariant parametrization",&
                       "Fock matrix parametrization",&
                       "Original matrix exponential parametrization"),&
         default_i_val=pao_rotinv_param)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    ! Preopt  ******************************************************************
    CALL keyword_create(keyword, name="EPS_PREOPT",&
         description="Convergence criteria for pre-optimization.",&
         default_r_val=1.e-6_dp)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="PREOPT_DM_FILE",&
         description="Read density matrix for preopt from given file.",&
         repeats=.FALSE., default_c_val="")
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    ! Misc ********************************************************************
    CALL keyword_create(keyword, name="READ_RESTART",&
         description="Reads given files as restart for PAO basis",&
         repeats=.FALSE., default_c_val="restart.pao")
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="CHECK_GRADIENT_PARAM_TOL",&
         description="Tolerance for check of analytic gradient of parametrization against the numeric one."//&
         " Negative values mean don't check at all.",&
         default_r_val=-1.0_dp)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="CHECK_GRADIENT_MCWEENY_TOL",&
         description="Tolerance for check of analytic gradient of McWeeny against the numeric one."//&
         " Negative values mean don't check at all.",&
         default_r_val=-1.0_dp)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="CHECK_UNITARY_TOL",&
         description="Check if rotation matrix is unitary."//&
         " Negative values mean don't check at all.",&
         default_r_val=-1.0_dp)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    ! Output *******************************************************************
    CALL section_create(subsection,name="PRINT",&
         description="Printkey section",&
         n_keywords=0, n_subsections=1, repeats=.TRUE.)

    CALL cp_print_key_section_create(printkey,"RUN_INFO",&
         description="Normal output by PAO", &
         print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__")
    CALL section_add_subsection(subsection,printkey)
    CALL section_release(printkey)

    CALL cp_print_key_section_create(printkey,"CG_INFO",&
         description="Output by the CG optimizer", &
         print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__")
    CALL section_add_subsection(subsection,printkey)
    CALL section_release(printkey)

    CALL cp_print_key_section_create(printkey,"RESTART",&
         description="Restart file of PAO basis", &
         print_level=high_print_level,add_last=add_last_numeric,filename="")
    CALL keyword_create(keyword, name="BACKUP_COPIES",&
         description="Specifies the maximum number of backup copies.",&
         usage="BACKUP_COPIES {int}",&
         default_i_val=1)
    CALL section_add_keyword(printkey,keyword)
    CALL keyword_release(keyword)

    CALL section_add_subsection(subsection,printkey)
    CALL section_release(printkey)

    CALL section_add_subsection(section,subsection)
    CALL section_release(subsection)

    ! CG stuff *****************************************************************
    CALL keyword_create(keyword, name="CG_INIT_STEPS",&
         description="Number of steepest descent steps before starting the"//&
         " conjugate gradients optimization.",&
         default_i_val=1)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL linesearch_create_section(subsection)
    CALL section_add_subsection(section, subsection)

    CALL section_release(subsection)
  END SUBROUTINE create_pao_section

END MODULE pao_input
