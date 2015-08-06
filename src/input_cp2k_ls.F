!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief input for the linear scaling (LS) section
!> \author Joost VandeVondele
! *****************************************************************************
MODULE input_cp2k_ls
  USE bibliography,                    ONLY: Lin2009,&
                                             Lin2013,&
                                             Niklasson2003,&
                                             Shao2003,&
                                             VandeVondele2012
  USE cp_output_handling,              ONLY: cp_print_key_section_create,&
                                             high_print_level
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE input_constants,                 ONLY: &
       ls_cluster_atomic, ls_cluster_molecular, ls_s_inversion_hotelling, &
       ls_s_inversion_sign_sqrt, ls_s_preconditioner_atomic, &
       ls_s_preconditioner_molecular, ls_s_preconditioner_none, &
       ls_scf_line_search_3point, ls_scf_line_search_3point_2d, ls_scf_ns, &
       ls_scf_pexsi, ls_scf_tc2, ls_scf_trs4
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: integer_t,&
                                             real_t
  USE kinds,                           ONLY: dp
  USE qs_density_mixing_types,         ONLY: create_mixing_section
  USE string_utilities,                ONLY: s2a
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_ls'

  PUBLIC :: create_ls_scf_section

CONTAINS
! *****************************************************************************
!> \brief creates the linear scaling scf section
!> \param section ...
!> \param error ...
!> \author Joost VandeVondele [2010-10]
! *****************************************************************************
  SUBROUTINE create_ls_scf_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ls_scf_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,"LS_SCF",&
         description="Specifies the parameters of the linear scaling SCF routines",&
         n_keywords=24, n_subsections=3, repeats=.FALSE., &
         citations=(/VandeVondele2012/),error=error)

    NULLIFY (keyword,subsection)

    CALL keyword_create(keyword, name="LS_DIIS",&
         description="Perform DIIS within linear scaling SCF",&
         usage="LS_DIIS",lone_keyword_l_val=.TRUE.,&
         default_l_val=.FALSE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="INI_DIIS",&
         description="Iteration cycle to start DIIS Kohn-Sham matrix update",&
         usage="INI_DIIS 2", default_i_val=2,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MAX_DIIS",&
         description="Size of LS_DIIS buffer",&
         usage="MAX_DIIS 4", default_i_val=4,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="NMIXING",&
         description="Minimal number of density mixing before start DIIS",&
         usage="NMIXING 2", default_i_val=2,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="EPS_DIIS",&
         description="Threshold on the convergence to start using DIIS",&
         usage="EPS_DIIS 1.e-1", default_r_val=1.e-1_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MAX_SCF",&
         description="Maximum number of SCF iteration to be performed for one optimization",&
         usage="MAX_SCF 200", default_i_val=20,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="EPS_SCF",&
         description="Target accuracy for the SCF convergence in terms of change of the total energy per electron.",&
         usage="EPS_SCF 1.e-6", default_r_val=1.e-7_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MIXING_FRACTION",&
         description="Mixing density matrices uses the specified fraction in the SCF procedure.",&
         usage="MIXING_FRACTION 0.4", default_r_val=0.45_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="EPS_FILTER",&
         description="Threshold used for filtering matrix operations.",&
         usage="EPS_FILTER 1.0E-7", default_r_val=1.0E-6_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="EPS_LANCZOS",&
         description="Threshold used for lanczos estimates.",&
         usage="EPS_LANCZOS 1.0E-4", default_r_val=1.0E-3_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MAX_ITER_LANCZOS",&
         description="Maximum number of lanczos iterations.",&
         usage="MAX_ITER_LANCZOS ", default_i_val=128,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MU",&
         description="Value (or initial guess) for the chemical potential,"//&
                     " i.e. some suitable energy between HOMO and LUMO energy.",&
         usage="MU 0.0", default_r_val=-0.1_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="FIXED_MU",&
         description="Should the calculation be performed at fixed chemical potential,"//&
                     " or should it be found fixing the number of electrons",&
         usage="FIXED_MU .TRUE.", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,  error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="EXTRAPOLATION_ORDER",&
         description="Number of previous matrices used for the ASPC extrapolation of the initial guess. "//&
                     "0 implies that an atomic guess is used at each step. "//&
                     "low (1-2) will result in a drift of the constant of motion during MD. "//&
                     "high (>5) might be somewhat unstable, leading to more SCF iterations.",&
         usage="EXTRAPOLATION_ORDER 3",default_i_val=4, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="S_PRECONDITIONER",&
         description="Preconditions S with some appropriate form.",&
         usage="S_PRECONDITIONER MOLECULAR",&
         default_i_val=ls_s_preconditioner_atomic,&
         enum_c_vals=s2a( "NONE", "ATOMIC", "MOLECULAR"),&
         enum_desc=s2a("No preconditioner",&
                       "Using atomic blocks",&
                       "Using molecular sub-blocks. Recommended if molecules are defined and not too large."),&
         enum_i_vals=(/ls_s_preconditioner_none, ls_s_preconditioner_atomic, ls_s_preconditioner_molecular/),&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="PURIFICATION_METHOD",&
         description="Scheme used to purify the Kohn-Sham matrix into the density matrix.",&
         usage="PURIFICATION_METHOD TRS4",&
         default_i_val=ls_scf_ns,&
         citations=(/VandeVondele2012,Niklasson2003/),&
         enum_c_vals=s2a( "SIGN_MATRIX", "TRS4", "TC2", "PEXSI" ),&
         enum_desc=s2a("Newton-Schulz sign matrix iteration.",&
                       "Trace resetting 4th order scheme",&
                       "Trace conserving 2nd order scheme",&
                       "PEXSI method"),&
         enum_i_vals=(/ls_scf_ns, ls_scf_trs4, ls_scf_tc2, ls_scf_pexsi/),&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="DYNAMIC_THRESHOLD",&
         description="Should the threshold for the purification be chosen dynamically",&
         usage="DYNAMIC_THRESHOLD .TRUE.", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,  error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="NON_MONOTONIC",&
         description="Should the purification be performed non-monotonically. Relevant for TC2 only.",&
         usage="NON_MONOTONIC .TRUE.", default_l_val=.TRUE., lone_keyword_l_val=.TRUE.,  error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MATRIX_CLUSTER_TYPE",&
         description="Specify how atomic blocks should be clustered in the used matrices, in order to improve flop rate, "//&
                     "and possibly speedup the matrix multiply. Note that the atomic s_preconditioner can not be used."//&
                     "Furthermore, since screening is on matrix blocks, "//&
                     "slightly more accurate results can be expected with molecular.",&
         usage="MATRIX_CLUSTER_TYPE MOLECULAR",&
         default_i_val=ls_cluster_atomic,&
         enum_c_vals=s2a("ATOMIC", "MOLECULAR"),&
         enum_desc=s2a("Using atomic blocks",&
                       "Using molecular blocks."),&
         enum_i_vals=(/ls_cluster_atomic, ls_cluster_molecular /),&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="SINGLE_PRECISION_MATRICES",&
         description="Matrices used within the LS code can be either double or single precision.",&
         usage="SINGLE_PRECISION_MATRICES", default_l_val=.FALSE.,lone_keyword_l_val=.TRUE., error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="RESTART_WRITE",&
         description="Write the density matrix at the end of the SCF (currently requires EXTRAPOLATION_ORDER>0). "//&
                     "Files might be rather large.",&
         usage="RESTART_READ", default_l_val=.FALSE.,lone_keyword_l_val=.TRUE., error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="RESTART_READ",&
         description="Read the density matrix before the (first) SCF.",&
         usage="RESTART_READ", default_l_val=.FALSE.,lone_keyword_l_val=.TRUE., error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="S_INVERSION",&
         description="Method used to compute the inverse of S.",&
         usage="S_PRECONDITIONER MOLECULAR",&
         default_i_val=ls_s_inversion_sign_sqrt,&
         enum_c_vals=s2a("SIGN_SQRT","HOTELLING"),&
         enum_desc=s2a("Using the inverse sqrt as obtained from sign function iterations.",&
                       "Using the Hotellign iteration."),&
         enum_i_vals=(/ls_s_inversion_sign_sqrt,ls_s_inversion_hotelling/),&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="SIGN_SQRT_ORDER",&
         description="Order of the sqrt iteration, should be 2..5, 3 or 5 recommended",&
         usage="SIGN_SQRT_ORDER 5",&
         default_i_val=3,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="REPORT_ALL_SPARSITIES",&
         description="Run the sparsity report at the end of the SCF",&
         usage="REPORT_ALL_SPARSITIES", default_l_val=.TRUE.,lone_keyword_l_val=.TRUE., error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="PERFORM_MU_SCAN",&
         description="Do a scan of the chemical potential after the SCF",&
         usage="PERFORM_MU_SCAN", default_l_val=.FALSE.,lone_keyword_l_val=.TRUE., error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL create_ls_curvy_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_chebyshev_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_mixing_section(subsection,ls_scf=.TRUE.,error=error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

    CALL create_pexsi_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)


  END SUBROUTINE create_ls_scf_section

! *****************************************************************************
!> \brief creates the DOS section
!> \param section ...
!> \param error ...
!> \author Joost VandeVondele, Jinwoong Cha [2012-10]
! *****************************************************************************
  SUBROUTINE create_chebyshev_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_chebyshev_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)

    CALL section_create(section,"CHEBYSHEV",&
         description="Specifies the parameters needed for the chebyshev expansion based properties.",&
         n_keywords=24, n_subsections=3, repeats=.FALSE., &
         error=error)

    NULLIFY (keyword)
    NULLIFY(print_key)

    CALL keyword_create(keyword, name="N_CHEBYSHEV",&
         description="Order of the polynomial expansion.",&
         usage="N_CHEBYSHEV 2000", default_i_val=500,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    ! A DOS print key
    CALL cp_print_key_section_create(print_key,"DOS",&
         description="Controls the printing of the Density of States (DOS).", &
         print_level=high_print_level,filename="",error=error)
    CALL keyword_create(keyword, name="N_GRIDPOINTS",&
         description="Number of points in the computed DOS",&
         usage="N_GRIDPOINTS 10000", default_i_val=2000,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL section_add_subsection(section, print_key, error=error)
    CALL section_release(print_key, error=error)

    ! Energy specific electron density cubes 
    CALL cp_print_key_section_create(print_key,&
         name="PRINT_SPECIFIC_E_DENSITY_CUBE",&
         description="Controls the printing of cube files with "//&
                     "the electronic density (states) "//&
                     "contributing to the density of states within "//&
                     "the specific energy range "//&
                     "(MIN_ENERGY &le; E &le; MAX_ENERGY). MIN_ENERGY and MAX_ENERGY need to be specified explicitly.",&
         print_level=high_print_level,filename="",error=error)

    CALL keyword_create(keyword, name="stride",&
         description="The stride (X,Y,Z) used to write the cube file "//&
         "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"//&
         " 1 number valid for all components.",&
         usage="STRIDE 2 2 2",n_var=-1,default_i_vals=(/2,2,2/), type_of_var=integer_t,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MIN_ENERGY",&
         description="Lower bounds of the energy ranges of interest.",&
         usage="MIN_ENERGY -1.01 -0.62 0.10 .. ", &
         type_of_var=real_t, n_var=-1, error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MAX_ENERGY",&
         description="Upper bounds of the energy ranges of interest.",&
         usage="MAX_ENERGY -0.81 -0.43 0.22 .. ", &
         type_of_var=real_t, n_var=-1, error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL section_add_subsection(section, print_key, error=error)
    CALL section_release(print_key, error=error)


  END SUBROUTINE create_chebyshev_section

! *****************************************************************************
!> \brief creates the curvy_steps section in linear scaling scf 
!> \param section ...
!> \param error ...
!> \author Florian Schiffmann [2012-10]
! *****************************************************************************
  SUBROUTINE create_ls_curvy_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ls_curvy_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,"CURVY_STEPS",&
         description="Specifies the parameters of the linear scaling SCF routines",&
         n_keywords=24, n_subsections=3, repeats=.FALSE., &
         citations=(/Shao2003/),&
         error=error)

    NULLIFY (keyword)

    CALL keyword_create(keyword, name="LINE_SEARCH",&
         description="Line serch type used in the curvy_setp optimization.",&
         usage="LINE Search 3POINT", default_i_val=ls_scf_line_search_3point,&
         enum_c_vals=s2a("3POINT","3POINT_2D"),&
         enum_desc=s2a("Performs a three point line search",&
         "Only for spin unrestricted calcualtions. Separate step sizes for alpha and beta spin"//&
         " using a fit to a 2D parabolic function"),&
         enum_i_vals=(/ls_scf_line_search_3point,ls_scf_line_search_3point_2d/), error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="N_BCH_HISTORY",&
         description="Number of stored matrices in the Baker-Campbell-Hausdorff series. "//&
         "Reduces the BCH evaluation during line search but can be memory intense. ",&
         usage="N_BCH_HISTORY 5",&
         default_i_val=7,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MIN_HESSIAN_SHIFT",&
         description="Minimal eigenvalue shift for the Hessian in the Newton iteration."//&
                     " Useful for small band gap systems (0.5-1.0 recommended). ",&
         usage="MIN_HESSIAN_SHIFT 0.0", default_r_val=0.0_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="FILTER_FACTOR",&
         description="Allows to set a seperate EPS_FILTER in the newton iterations."//&
                     " The new EPS is EPS_FILTER*FILTER_FACTOR.",&
         usage="FILTER_FACTOR 10.0", default_r_val=1.0_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="FILTER_FACTOR_SCALE",&
         description="Allows for dynamic EPS_FILTER. Updates the filter factor every SCF-Newton "//&
                     "step by FILTER_FACTOR=FILTER_FACTOR*FILTER_FACTOR_SCALE",&
         usage="FILTER_FACTOR_SCALE 0.5", default_r_val=1.0_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MIN_FILTER",&
         description="Lowest EPS_FILTER in dynamic filtering. Given as multiple of EPS_FILTER:"//&
                     " EPS_FILTER_MIN=EPS_FILTER*MIN_FILTER",&
         usage="FILTER_FACTOR 1.0", default_r_val=1.0_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)


  END SUBROUTINE create_ls_curvy_section

! *****************************************************************************
!> \brief creates the PEXSI library subsection of the linear scaling section.
!> \param section ...
!> \param error ...
!> \par History
!>      11.2014 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  SUBROUTINE create_pexsi_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_pexsi_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)

    CALL section_create(section,"PEXSI",&
         description="Specifies the parameters of the PEXSI library. The density matrix is calculated "//&
         "with PEXSI if PURIFICATION_METHOD (in LS_SCF section) is set to PEXSI. "//&
         "The computational cost of PEXSI is at most quadratically scaling w.r.t. the system size "//& 
         "and PEXSI is applicable to insulating and metallic systems. The value of EPS_PGF_ORB "//&
         "(in QS input section) defines the sparsity of the matrices sent to PEXSI and EPS_FILTER "//&
         "is overwritten with 0.",&
         n_keywords=17, repeats=.FALSE.,citations=(/Lin2009,Lin2013/),error=error)
     NULLIFY (keyword)

     CALL keyword_create(keyword, name="TEMPERATURE",&
          description="Electronic temperature",&
          default_r_val=cp_unit_to_cp2k(value=300.0_dp, unit_str="K", error=error),&
          unit_str="K", error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="GAP",&
          description="Spectral gap. Note: This can be set to be 0 in most cases.",&
          default_r_val=0.0_dp,unit_str="hartree",error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="NUM_POLE",&
          description="Number of terms in the pole expansion (should be even).",&
          default_i_val=64,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="IS_INERTIA_COUNT",&
          description="Whether inertia counting is used each time the DFT driver "//& 
          "of PEXSI is invoked. If FALSE, inertia counting is still used in the "//&
          "first SCF iteration.",&
          default_l_val=.FALSE., lone_keyword_l_val = .TRUE. ,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="MAX_PEXSI_ITER",&
          description="Maximum number of PEXSI iterations after each inertia counting procedure.",&
          default_i_val=5,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="MU_MIN_0",&
          description="Initial guess of lower bound for mu.",&
          default_r_val=-5.0_dp,unit_str="hartree",error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="MU_MAX_0",&
          description="Initial guess of upper bound for mu.",&
          default_r_val=5.0_dp,unit_str="hartree",error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="MU_INERTIA_TOLERANCE",&
          description="Stopping criterion in terms of the chemical potential for the "//&
          "inertia counting procedure.",&
          default_r_val=0.01_dp,unit_str="hartree",error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="MU_INERTIA_EXPANSION",&
          description="If the chemical potential is not in the initial interval, "//&
          "the interval is expanded by MU_INERTIA_EXPANSION.",&
          default_r_val=0.15_dp,unit_str="hartree",error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="MU_PEXSI_SAFE_GUARD",&
          description="Safe guard criterion in terms of the chemical potential to "//& 
          "reinvoke the inertia counting procedure.",&
          default_r_val=0.01_dp,unit_str="hartree",error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="NUM_ELECTRON_PEXSI_TOLERANCE",&
          description="Stopping criterion of the PEXSI iteration in terms of "//&
          "The number of electrons compared to the exact number of electrons. "//&
          "This threshold is the target tolerance applied at convergence of SCF.",&
          default_r_val=0.1_dp,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="NUM_ELECTRON_INITIAL_TOLERANCE",&
          description="The same as NUM_ELECTRON_PEXSI_TOLERANCE but applied in the first SCF steps. "//&
                      "If set to a value smaller than NUM_ELECTRON_PEXSI_TOLERANCE, it is overwritten with "//&
                      "NUM_ELECTRON_PEXSI_TOLERANCE (default). If set to a value larger than "//&
                      "NUM_ELECTRON_PEXSI_TOLERANCE, the PEXSI tolerance in number of electrons is set adaptively "//&
                      "according to the SCF convergence error of the previous SCF step. This reduces the number "//&
                      "of PEXSI iterations in the first SCF steps but leads to at least one more SCF step.",&
          default_r_val=0.0_dp,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="ORDERING",&
          description="Ordering strategy for factorization and selected inversion.",&
          enum_c_vals=s2a("PARALLEL","SEQUENTIAL","MULTIPLE_MINIMUM_DEGREE"),&
          enum_desc=s2a("Parallel ordering using ParMETIS/PT-SCOTCH (PARMETIS option in SuperLU_DIST)",&
          "Sequential ordering using METIS (METIS_AT_PLUS_A option in SuperLU_DIST)",&
          "Multiple minimum degree ordering (MMD_AT_PLUS_A option in SuperLU_DIST)"),&
          enum_i_vals=(/0,1,2/), default_i_val=0, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="NP_SYMB_FACT",&
          description="Number of processors for PARMETIS/PT-SCOTCH. Only used if the ordering == 0. "//&
          "If 0, the number of processors for PARMETIS/PT-SCOTCH will be set equal to the number of "//&
          "MPI ranks per pole. Note: if more than one processor is used, a segmentation fault may occur in the "//&
          "symbolic factorization phase.",&
          default_i_val=1,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="VERBOSITY",&
          description="The level of output information.",&
          enum_c_vals=s2a("SILENT","BASIC","DETAILED"),&
          enum_i_vals=(/0,1,2/), default_i_val=1, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="MIN_RANKS_PER_POLE",&
          description="The minimum number of processors used for each pole. The real "//&
          "number of processors per pole is the smallest number greater or equal to "//&
          "MIN_RANKS_PER_POLE that divides MPI size without remainder. For efficiency, MIN_RANKS_PER_POLE "//&
          "should be a small numbers (limited by the available memory).", &
          default_i_val=64,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="CSR_SCREENING",&
          description="Whether distance screening should be applied to improve sparsity of CSR matrices.",&
          default_l_val=.TRUE., lone_keyword_l_val = .TRUE. ,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

  END SUBROUTINE create_pexsi_section

END MODULE input_cp2k_ls
