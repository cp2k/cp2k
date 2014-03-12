!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief input for the ALMO SCF section
!> \author Rustam Khaliullin
! *****************************************************************************
MODULE input_cp2k_almo
  USE bibliography,                    ONLY: Khaliullin2013
  USE f77_blas
  USE input_constants
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_create,&
                                             section_type
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_almo'

  PUBLIC :: create_almo_scf_section

CONTAINS

! *****************************************************************************
!> \brief create the almo scf section
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
SUBROUTINE create_almo_scf_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_almo_scf_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"ALMO_SCF",&
            description="Settings for methods based on absolutely localized orbitals"//&
            " (ALMOs). Restricted to closed-shell molecular systems.",&
            n_keywords=11, n_subsections=0, repeats=.FALSE.,&
            citations=(/Khaliullin2013/),&
            required=.FALSE., error=error)

       NULLIFY (keyword)

       CALL keyword_create(keyword, name="BLOCKED_MAX_ITER",&
            description="Maximum number of SCF iterations for optimization of block-diagonal ALMOs",&
            usage="BLOCKED_MAX_ITER 200", default_i_val=100,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BLOCKED_N_DIIS",&
            description="Number of error vectors to be used in the DIIS "//&
            "optimization procedure of block-diagonal ALMOs",&
            usage="BLOCKED_N_DIIS 5", default_i_val=7,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BLOCKED_EPS_ITER",&
            description="Target accuracy for the convergence of optimization of block-diagonal ALMOs",&
            usage="BLOCKED_EPS_ITER 1.e-5", default_r_val=1.e-5_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DELOCALIZE_MAX_ITER",&
            description="Maximum number of SCF iterations for optimization of ALMOs",&
            usage="DELOCALIZE_MAX_ITER 200", default_i_val=100,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DELOCALIZE_EPS_ITER",&
            description="Target accuracy for the convergence of optimization of ALMOs",&
            usage="DELOCALIZE_EPS_ITER 1.e-5", default_r_val=1.e-5_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DELOCALIZE_EPS_LIN_SEARCH",&
            description="Target accuracy of the line search for the PCG optimization of ALMOs",&
            usage="DELOCALIZE_EPS_LIN_SEARCH 1.e-6", default_r_val=1.e-7_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_FILTER",&
            description="Threshold for the sparsity filter",&
            usage="EPS_FILTER 1.e-6", default_r_val=1.e-7_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="DELOCALIZE_CONJUGATOR",&
            description="Various methods to compute step directions in the PCG optimization of ALMOs",&
            usage="DELOCALIZE_CONJUGATOR POLAK_RIBIERE",&
            default_i_val=cg_hager_zhang,&
            enum_c_vals=s2a("ZERO", "POLAK_RIBIERE", "FLETCHER_REEVES",&
            "HESTENES_STIEFEL", "FLETCHER", "LIU_STOREY", "DAI_YUAN","HAGER_ZHANG"),&
            enum_desc=s2a("Steepest descent","Polak and Ribiere",&
            "Fletcher and Reeves","Hestenes and Stiefel",&
            "Fletcher (Conjugate descent)","Liu and Storey",&
            "Dai and Yuan","Hager and Zhang"),&
            enum_i_vals=(/cg_zero,cg_polak_ribiere,cg_fletcher_reeves,&
                          cg_hestenes_stiefel,cg_fletcher,cg_liu_storey,&
                          cg_dai_yuan,cg_hager_zhang/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DELOCALIZE_ALGORITHM",&
            description="Specifies the algorithm to update ALMOs.",&
            usage="DELOCALIZE_ALGORITHM BDIAG",&
            default_i_val=almo_scf_diag,&
            enum_c_vals=s2a("DDIAG", "PCG"),&
            enum_desc=s2a("Domain diagonalization",&
                          "Preconditioned conjugate gradient"),&
            enum_i_vals=(/almo_scf_diag,almo_scf_pcg/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DELOCALIZE_R_CUTOFF_FACTOR",&
            description="Parameter to calculate the localization radius: "//&
            !"r0 = r0_factor*(radius(at1)+radius(at2)) + r0_shift",&
            "R_cutoff = DELOCALIZE_R_CUTOFF_FACTOR*(radius(at1)+radius(at2))",&
            usage="DELOCALIZE_R_CUTOFF_FACTOR 1.6", default_r_val=1.60_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DELOCALIZE_METHOD",&
            description="Methods to treat electron delocalization after block-diagonal ALMOs are computed",&
            usage="DELOCALIZE_METHOD X_R",&
            default_i_val=almo_deloc_qx,&
            enum_c_vals=s2a("NONE","X","SCF","X_R","SCF_R","X_THEN_SCF"),&
            enum_desc=s2a("Neglect electron delocalization",&
            "Single excitation delocalization energy correction",&
            "Fully self-consistent treatment of delocalization",&
            "Spacially restricted single excitation delocalization energy correction",&
            "Spacially restricted self-consistent treatment of delocalization",&
            "Single excitation delocalization energy correction followed by full SCF procedure"),&
            enum_i_vals=(/almo_deloc_none,almo_deloc_x,almo_deloc_full_scf,&
            almo_deloc_qx,almo_deloc_qscf,&
            !almo_deloc_xk,&
            almo_deloc_x_then_full_scf/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    
    END IF

  END SUBROUTINE create_almo_scf_section

! *****************************************************************************
!> \brief create the optimizer section
!> \par History
!>       2012.03 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
SUBROUTINE create_optimizer_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_optimizer_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
    END IF

  END SUBROUTINE create_optimizer_section

END MODULE input_cp2k_almo
