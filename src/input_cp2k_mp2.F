!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief input section for MP2
!> \par History
!>      05.2011 created
!> \author MDB
! *****************************************************************************
MODULE input_cp2k_mp2
  USE bibliography,                    ONLY: DelBen2012,&
                                             DelBen2013,&
                                             DelBen2015
  USE cp_output_handling,              ONLY: add_last_numeric,&
                                             cp_print_key_section_create,&
                                             debug_print_level,&
                                             high_print_level,&
                                             low_print_level,&
                                             medium_print_level,&
                                             silent_print_level
  USE input_constants,                 ONLY: &
       do_hfx_potential_coulomb, do_mp2_potential_tshpsc, gaussian, &
       mp2_method_direct, mp2_method_gpw, mp2_method_none, &
       mp2_ri_optimize_basis, numerical, ri_mp2_laplace, ri_mp2_method_gpw, &
       ri_rpa_g0w0_crossing_none, ri_rpa_g0w0_crossing_z_shot, &
       ri_rpa_method_gpw, wfc_mm_style_gemm, wfc_mm_style_syrk
  USE input_cp2k_hfx,                  ONLY: create_hfx_section
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
  USE string_utilities,                ONLY: s2a
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_mp2'

  PUBLIC :: create_mp2_section

CONTAINS

! *****************************************************************************
!> \brief creates the input section for the mp2 part
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author MDB
! *****************************************************************************
  SUBROUTINE create_mp2_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mp2_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,"WF_CORRELATION",&
         description="Sets up the Wavefunction-based Correlation parameters if requested ",&
         n_keywords=8, n_subsections=6, repeats=.TRUE., &
         citations=(/DelBen2012,DelBen2013,DelBen2015/),&
         error=error)

    NULLIFY(keyword, print_key, subsection)

    CALL keyword_create(&
         keyword=keyword,&
         name="METHOD",&
         citations=(/DelBen2012/),&
         description="Which method should be used to compute the MP2 energy",&
         usage="METHOD MP2_GPW",&
         enum_c_vals=s2a("NONE","DIRECT_CANONICAL","MP2_GPW","RI_MP2_GPW","RI_RPA_GPW","RI_SOS_LAPLACE",&
                         "OPTIMIZE_RI_BASIS"),&
         enum_i_vals=(/mp2_method_none,mp2_method_direct,mp2_method_gpw,ri_mp2_method_gpw,&
                       ri_rpa_method_gpw,ri_mp2_laplace,mp2_ri_optimize_basis/),&
         enum_desc=s2a("Skip MP2 calculation",&
                       "Use the direct mp2 canonical approach",&
                       "Use the GPW approach to MP2",&
                       "Use the GPW approach to RI-MP2",&
                       "Use the GPW approach to RI-RPA",&
                       "Use the GPW approach to RI-SOS-Laplace-MP2",&
                       "Optimize RIMP2 basis set"),&
         default_i_val=mp2_method_direct, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(&
         keyword=keyword,&
         name="MEMORY",&
         description="Maximum allowed total memory usage during MP2 methods [MiB].",&
         usage="MEMORY 1500 ",&
         default_r_val=1.024E+3_dp,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(&
         keyword=keyword,&
         name="SCALE_S",&
         description="Scaling factor of the singlet energy component (opposite spin, OS). ",&
         usage="SCALE_S  1.0",&
         default_r_val=1.0_dp,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(&
         keyword=keyword,&
         name="SCALE_T",&
         description="Scaling factor of the triplet energy component (same spin, SS).",&
         usage="SCALE_T  1.0",&
         default_r_val=1.0_dp,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(&
         keyword=keyword,&
         name="GROUP_SIZE",&
         variants=(/"NUMBER_PROC"/),&
         description="Group size used in the computation of the integrals. Default is to use all processors (GROUP_SIZE=-1)."//&
          "A smaller group size (for example the node size), might a better choice if the actual MP2 time is large "//&
          "compared to integral computation time. This is usually the case if the total number of processors is not too large.",&
         usage="GROUP_SIZE 32",&
         default_i_val=-1,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(&
         keyword=keyword,&
         name="ROW_BLOCK",&
         variants=(/"ROW_BLOCK_SIZE"/),&
         description="Size of the row block used in the SCALAPACK block cyclic data distribution."//&
          "Default is (ROW_BLOCK=-1) is automatic. "//&
          "A proper choice can speedup the parallel matrix multiplication in the case of RI-RPA and RI-SOS-MP2-Laplace.",&
         usage="ROW_BLOCK 512",&
         default_i_val=-1,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(&
         keyword=keyword,&
         name="COL_BLOCK",&
         variants=(/"COL_BLOCK_SIZE"/),&
         description="Size of the column block used in the SCALAPACK block cyclic data distribution."//&
          "Default is (COL_BLOCK=-1) is automatic. "//&
          "A proper choice can speedup the parallel matrix multiplication in the case of RI-RPA and RI-SOS-MP2-Laplace.",&
         usage="COL_BLOCK 512",&
         default_i_val=-1,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(&
         keyword=keyword,&
         name="CALC_COND_NUM",&
         variants=(/"CALC_CONDITION_NUMBER"/),&
         description="Calculate the condition number of the (P|Q) matrix for the RI methods.",&
         usage="CALC_COND_NUM",&
         default_l_val=.FALSE.,&
         lone_keyword_l_val=.TRUE.,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL cp_print_key_section_create(print_key,"MP2_INFO",&
         description="Controls the printing basic info about MP2 method", &
         print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
         error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL create_mp2_direct(subsection,error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

    CALL create_wfc_gpw(subsection,error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

    CALL create_ri_mp2(subsection,error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

    CALL create_opt_ri_basis(subsection,error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

    CALL create_ri_rpa(subsection,error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

    CALL create_ri_laplace(subsection,error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

    CALL create_cphf(subsection,error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

    CALL create_mp2_potential(subsection,error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)


  END SUBROUTINE create_mp2_section

! *****************************************************************************
!> \brief ...
!> \param section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_mp2_direct(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mp2_direct', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,"DIRECT_CANONICAL",&
         description="Parameters influencing the direct canonical method",&
         n_keywords=1, n_subsections=0, repeats=.FALSE., &
         error=error)

    NULLIFY(keyword)

    CALL keyword_create(&
         keyword=keyword,&
         name="BIG_SEND",&
         description="Send big messages between processes (useful for >48 processors).",&
         usage="BIG_SEND",&
         default_l_val=.TRUE.,&
         lone_keyword_l_val=.TRUE.,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)


  END SUBROUTINE create_mp2_direct

! *****************************************************************************
!> \brief ...
!> \param section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_ri_mp2(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ri_mp2', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,"RI_MP2",&
         description="Parameters influencing the RI MP2 method",&
         n_keywords=3, n_subsections=0, repeats=.FALSE., &
         citations=(/DelBen2013/),&
         error=error)

    NULLIFY(keyword)

    CALL keyword_create(keyword, name="BLOCK_SIZE",&
         variants=(/"MESSAGE_SIZE"/),&
         description="Determines the blocking used for communication in RI-MP2. Larger BLOCK_SIZE "//&
                     "reduces communication but requires more memory. The default (-1) is automatic.",&
         usage="BLOCK_SIZE 2",&
         default_i_val=-1,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="EPS_CANONICAL",&
         description="Threshold for discriminate if a given ij pairs "//&
                     "of the unrelaxed MP2 density matrix has to be "//&
                     "calculated with a canonical reformulation based "//&
                     "on the occupied eigenvalues differences.",&
         usage="EPS_CANONICAL 1.0E-8",type_of_var=real_t,&
         default_r_val=1.0E-7_dp, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(&
         keyword=keyword,&
         name="FREE_HFX_BUFFER",&
         description="Free the buffer containing the 4 center integrals used in the Hartree-Fock exchange calculation. "//&
                     "This will be effective only for gradients calculations, since for the energy only "//&
                     "case, the buffers are released by default. (Right now debugging only).",&
         usage="FREE_HFX_BUFFER",&
         default_l_val=.TRUE.,&
         lone_keyword_l_val=.TRUE.,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)


  END SUBROUTINE create_ri_mp2

! *****************************************************************************
!> \brief ...
!> \param section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_opt_ri_basis(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_opt_ri_basis', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,"OPT_RI_BASIS",&
         description="Parameters influencing the optimization of the RI MP2 basis. "//&
                     "Only exponents of non-contracted auxiliary basis can be optimized. "//&
                     "An initial RI auxiliary basis has to be specified.",&
         n_keywords=6, n_subsections=0, repeats=.FALSE., &
         citations=(/DelBen2013/),&
         error=error)

    NULLIFY(keyword)

    CALL keyword_create(keyword, name="DELTA_I_REL",&
         variants=(/"DI_REL"/),&
         description="Target accuracy in the relative deviation of the amplitudes calculated with "//&
                     "and without RI approximation, (more details in Chem.Phys.Lett.294(1998)143).",&
         usage="DELTA_I_REL  1.0E-6_dp",&
         default_r_val=1.0E-6_dp,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="DELTA_RI",&
         variants=(/"DRI"/),&
         description="Target accuracy in the absolute difference between the RI-MP2 "//&
                     "and the exact MP2 energy, DRI=ABS(E_MP2-E_RI-MP2).",&
         usage="DELTA_RI  1.0E-6_dp",&
         default_r_val=5.0E-6_dp,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="EPS_DERIV",&
         variants=(/"EPS_NUM_DERIV"/),&
         description="The derivatives of the MP2 energy with respect to the "//&
                     "exponents of the basis are calculated numerically. "//&
                     "The change in the exponent a_i employed for the numerical evaluation "//&
                     "is defined as h_i=EPS_DERIV*a_i.",&
         usage="EPS_DERIV  1.0E-3_dp",&
         default_r_val=1.0E-3_dp,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MAX_ITER",&
         variants=(/"MAX_NUM_ITER"/),&
         description="Specifies the maximum number of steps in the RI basis optimization.",&
         usage="MAX_ITER 100",&
         default_i_val=50,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="NUM_FUNC",&
         description="Specifies the number of function, for each angular momentum (s, p, d ...), "//&
                     "employed in the automatically generated initial guess. "//&
                     "This will be effective only if RI_AUX_BASIS_SET in the KIND section is not specified.",&
         usage="NUM_FUNC {number of s func.} {number of p func.} ...", &
         n_var=-1, default_i_vals=(/-1/), type_of_var=integer_t, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword=keyword,name="BASIS_SIZE",&
         description="Specifies the size of the auxiliary basis set automatically "//&
                     "generated as initial guess. This will be effective only if RI_AUX_BASIS_SET "//&
                     "in the KIND section and NUM_FUNC are not specified.",&
         usage="BASIS_SIZE  (MEDIUM|LARGE|VERY_LARGE)",&
         enum_c_vals=s2a("MEDIUM","LARGE","VERY_LARGE"),&
         enum_i_vals=(/0,1,2/),&
         default_i_val=0, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)


  END SUBROUTINE create_opt_ri_basis

! *****************************************************************************
!> \brief ...
!> \param section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_ri_laplace(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ri_laplace', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,"RI_LAPLACE",&
         description="Parameters influencing the RI-SOS-MP2-Laplace method",&
         n_keywords=2, n_subsections=0, repeats=.FALSE., &
         citations=(/DelBen2013/),&
         error=error)

    NULLIFY(keyword)

    CALL keyword_create(keyword, name="QUADRATURE_POINTS",&
         variants=(/"LAPLACE_NUM_QUAD_POINTS"/),&
         description="Number of quadrature points for the numerical integration in the RI-SOS-MP2-Laplace method.",&
         usage="QUADRATURE_POINTS 6",&
         default_i_val=5,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="SIZE_INTEG_GROUP",&
         variants=(/"LAPLACE_GROUP_SIZE"/),&
         description="Group size for the integration in the Laplace method, that is the number of processes involved in "//&
                     "the computation of each integration point. SIZE_INTEG_GROUP has to be a multiple "//&
                     "of GROUP_SIZE in the WF_CORRELATION section. The default (-1) "//&
                     "is automatic.",&
         usage="SIZE_INTEG_GROUP 16",&
         default_i_val=-1,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)


  END SUBROUTINE create_ri_laplace

! *****************************************************************************
!> \brief ...
!> \param section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_ri_rpa(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ri_rpa', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,"RI_RPA",&
         description="Parameters influencing the RI RPA method",&
         n_keywords=2, n_subsections=1, repeats=.FALSE., &
         citations=(/DelBen2013,DelBen2015/),&
         error=error)

    NULLIFY(keyword,subsection)

    CALL keyword_create(keyword, name="QUADRATURE_POINTS",&
         variants=(/"RPA_NUM_QUAD_POINTS"/),&
         description="Number of quadrature points for the numerical integration in the RI-RPA method.",&
         usage="QUADRATURE_POINTS 60",&
         default_i_val=40,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="SIZE_FREQ_INTEG_GROUP",&
         variants=(/"RPA_GROUP_SIZE"/),&
         description="Group size for frequency integration, that is the number of processes involved in "//&
                     "the computation of each integration point. SIZE_FREQ_INTEG_GROUP has to be a multiple "//&
                     "of GROUP_SIZE in the WF_CORRELATION section. The default (-1) "//&
                     "is automatic.",&
         usage="SIZE_FREQ_INTEG_GROUP 16",&
         default_i_val=-1,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword,&
         name="MM_STYLE",&
         description="Matrix multiplication style for the Q matrix.",&
         usage="MM_STYLE GEMM",&
         enum_c_vals=s2a("GEMM","SYRK"),&
         enum_i_vals=(/wfc_mm_style_gemm,wfc_mm_style_syrk/),&
         enum_desc=s2a("Use pdgemm: more flops, maybe faster.",&
                       "Use pdysrk: fewer flops, maybe slower."),&
         default_i_val=wfc_mm_style_gemm, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(&
         keyword=keyword,&
         name="MINIMAX_QUADRATURE",&
         variants=(/"MINIMAX"/),&
         description="Use the Minimax quadrature scheme for performing the numerical integration. "//&
                     "Maximum number of quadrature point limited to 20.",&
         usage="MINIMAX_QUADRATURE",&
         default_l_val=.FALSE.,&
         lone_keyword_l_val=.TRUE.,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(&
         keyword=keyword,&
         name="RI_G0W0",&
         variants=(/"GW"/),&
         description="Decide whether to perform an RI_G0W0 calculation on top of RI_ RPA.",&
         usage="RI_G0W0",&
         default_l_val=.FALSE.,&
         lone_keyword_l_val=.TRUE.,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)


    ! here we generate a hfx subsection to use in the case EXX has to be computed after RPA
    CALL create_hfx_section(subsection,error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

    ! here we generate a G0W0 subsection to use if G0W0 is desired
    CALL create_ri_g0w0(subsection,error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)


  END SUBROUTINE create_ri_rpa

! *****************************************************************************
!> \brief ...
!> \param section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_ri_g0w0(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ri_g0w0', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,"RI_G0W0",&
         description="Parameters influencing the RI-G0W0 method",&
         n_keywords=12, n_subsections=0, repeats=.FALSE., &
         citations=(/DelBen2013/),&
         error=error)

    NULLIFY(keyword)

    CALL keyword_create(keyword, name="CORR_MOS_OCC",&
         variants=(/"CORR_OCC"/),&
         description="Number of occupied MOs whose energies are corrected by RI-G0W0. "//&
                     "Counting beginning from HOMO, e.g. 3 corrected occ. MOs correspond "//&
                     "to correction of HOMO, HOMO-1 and HOMO-2. Numerical effort and "//&
                     "storage of RI-G0W0 increase linearly with this number. In case you "//&
                     "want to correct all occ. MOs, insert a number larger than the number "//&
                     "of occ. MOs.",&
         usage="CORR_OCC 3",&
         default_i_val=10,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="CORR_MOS_VIRT",&
         variants=(/"CORR_VIRT"/),&
         description="Number of virtual MOs whose energies are corrected by RI-G0W0. "//&
                     "Counting beginning from LUMO, e.g. 3 corrected occ. MOs correspond "//&
                     "to correction of LUMO, LUMO+1 and LUMO+2. Numerical effort and "//&
                     "storage of RI-G0W0 increase linearly with this number. In case you "//&
                     "want to correct all virt. MOs, insert a number larger than the number "//& 
                     "of virt. MOs.",&
         usage="CORR_VIRT 3",&
         default_i_val=10,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="SCALING",&
         variants=(/"A_SCALING"/),&
         description="Set the scaling factor a for the Clenshaw-Curtis grid by hand. "//&
                     "If SCALING < 0, a is calculated as for RI-RPA. Until now, the best "//&
                     "choice seems to be 0.2. The grid points of the frequency grid are "//&
                     "just multiplied by SCALING: omega_j = omega_j*Scaling.",&
         usage="SCALING 0.2",&
         default_r_val=0.2_dp,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="NUMB_POLES",&
         description="Number of poles for the fitting. Usually, two poles are sufficient. ",&
         usage="NUMB_POLES 2",&
         default_i_val=2,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="OMEGA_MAX_FIT",&
         description="Determines fitting range for the self-energy on the imaginary axis: "//&
                     "[0, OMEGA_MAX_FIT] for virt orbitals, [-OMEGA_MAX_FIT,0] for occ orbitals. "//&
                     "Unit: Hartree. Default: 0.734996 H = 20 eV. ",&
         usage="OMEGA_MAX_FIT 0.5",&
         default_r_val=0.734996_dp,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="STOP_CRIT",&
         variants=(/"STOP_CRIT_1"/),&
         description="Convergence criterion for the fit. The fitting iteration terminates, if "//&
                     "chi^2(n+1)/chi^2(n)-1 < STOP_CRIT.",&
         usage="STOP_CRIT 1.0E-7",&
         default_r_val=1E-5_dp,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="PRINT_FIT_ERROR",&
         variants=(/"FIT_ERROR"/),&
         description="Print the statistical error due to the fitting. This error in most cases "//&
                     "is very pessimistic (e.g. can be a factor of 10 too large).",&
         usage="FIT_ERROR",&
         default_l_val=.FALSE.,&
         lone_keyword_l_val=.TRUE.,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MAX_ITER_FIT",&
         description="Maximum number of iterations for the Levenberg-Marquard fit.",&
         usage="MAX_ITER_FIT 10000",&
         default_i_val=10000,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="CHECK_FIT",&
         description="If true, the self-energy on the imaginary axis and the fit are printed"//&
                     "the file self_energy_of_MO_<level>_for_imaginary_frequency.",&
         usage="CHECK_FIT",&
         default_l_val=.FALSE.,&
         lone_keyword_l_val=.TRUE.,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="CROSSING_SEARCH",&
         description="Determines, how the self_energy is evaluated on the real axis.",&
         usage="CROSSING_SEARCH Z_SHOT",&
         enum_c_vals=s2a("NONE","Z_SHOT"),&
         enum_i_vals=(/ri_rpa_g0w0_crossing_none,ri_rpa_g0w0_crossing_z_shot/),&
         enum_desc=s2a("Evaluate the correlation self-energy at the energy eigenvalue of SCF.",&
                       "Calculate the derivative of Sigma and out of it Z. Then extrapolate using Z"),&
         default_i_val=ri_rpa_g0w0_crossing_z_shot, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="FERMI_LEVEL_OFFSET",&
          description="Fermi level for occ. orbitals: e_HOMO + FERMI_LEVEL_OFFSET; "//&
                      "Fermi level for virt. orbitals: e_LUMO - FERMI_LEVEL_OFFSET. "//&
                      "In case e_homo + FERMI_LEVEL_OFFSET < e_lumo - FERMI_LEVEL_OFFSET,"//&
                      "we set Fermi level = (e_HOMO+e_LUMO)/2",&
          usage="FERMI_LEVEL_OFFSET 1.0E-2",&
          default_r_val=2.0E-2_dp,&
          error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="CUTOFF_RADIUS",&
          description="Cutoff radius for the truncated Coulomb potential in GW in Angstrom. "//&
                      "The cell length seems to be a good choice for the truncation.",&
          usage="CUTOFF_RADIUS 6.0",&
          default_r_val=10.0_dp,&
          unit_str="angstrom",error=error)


    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="TRUNCATION",&
         description="If true, the Coulomb potential is truncated. Recommended for periodic"//&
                     "systems.",&
         usage="TRUNCATION",&
         default_l_val=.FALSE.,&
         lone_keyword_l_val=.TRUE.,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="EV_SC_ITER",&
         description="Number of iterations for eigenvalue self-consistency cycle. The "//&
                     "computational effort of RPA and GW scales linearly with this number.",&
         usage="EV_SC_ITER 3",&
         default_i_val=1,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="HF_LIKE_EV_START",&
         description="If true, take as input for GW/RPA corrected HF-like eigenvalues according "//&
                     "to PRB 83, 115103 (2011), Sec. IV.",&
         usage="HF_LIKE_EV_START",&
         default_l_val=.FALSE.,&
         lone_keyword_l_val=.TRUE.,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="PRINT_GW_DETAILS",&
         description="If true, prints additional information on the quasiparticle energies.",&
         usage="PRINT_GW_DETAILS",&
         default_l_val=.FALSE.,&
         lone_keyword_l_val=.TRUE.,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)


  END SUBROUTINE create_ri_g0w0

! *****************************************************************************
!> \brief ...
!> \param section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_wfc_gpw(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_wfc_gpw', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,"WFC_GPW",&
         description="Parameters for the GPW approach in Wavefunction-based Correlation methods",&
         n_keywords=5, n_subsections=0, repeats=.FALSE., &
         error=error)

    NULLIFY(keyword)
    CALL keyword_create(keyword, name="EPS_GRID",&
         description="Determines a threshold for the GPW based integration",&
         usage="EPS_GRID 1.0E-9 ",type_of_var=real_t,&
         default_r_val=1.0E-8_dp, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="EPS_FILTER",&
         description="Determines a threshold for the DBCSR based multiply (usually 10 times smaller than EPS_GRID).",&
         usage="EPS_FILTER 1.0E-10 ",type_of_var=real_t,&
         default_r_val=1.0E-9_dp, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="CUTOFF",&
         description="The cutoff of the finest grid level in the MP2 gpw integration.",&
         usage="CUTOFF 300",type_of_var=real_t,&
         default_r_val=300.0_dp, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="REL_CUTOFF",&
         variants=(/"RELATIVE_CUTOFF"/),&
         description="Determines the grid at which a Gaussian is mapped.",&
         usage="REL_CUTOFF 50",type_of_var=real_t,&
         default_r_val=50.0_dp, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MULTIPOLE_TWO_CENT_INT",&
         description="Use fast algorithm to calculate two-center RI integrals based on multipole expansion",&
         usage="MULTIPOLE_TWO_CENT_INT",&
         default_l_val=.FALSE.,&
         lone_keyword_l_val=.TRUE.,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="PRINT_LEVEL",&
         variants=(/"IOLEVEL"/),&
         description="How much output is written by the individual groups.",&
         usage="PRINT_LEVEL HIGH",&
         default_i_val=silent_print_level,enum_c_vals=&
          s2a("SILENT","LOW","MEDIUM","HIGH","DEBUG"),&
         enum_desc=s2a( "Almost no output",&
                        "Little output", "Quite some output", "Lots of output",&
                        "Everything is written out, useful for debugging purposes only"),&
         enum_i_vals=(/silent_print_level,low_print_level,medium_print_level,&
         high_print_level,debug_print_level/),error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)


  END SUBROUTINE create_wfc_gpw

! *****************************************************************************
!> \brief ...
!> \param section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_cphf(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_cphf', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,"CPHF",&
         description="Parameters influencing the solution of the Z-vector equations in MP2 gradients calculations.",&
         n_keywords=2, n_subsections=0, repeats=.FALSE., &
         citations=(/DelBen2013/),&
         error=error)

    NULLIFY(keyword)

    CALL keyword_create(keyword, name="MAX_ITER",&
         variants=(/"MAX_NUM_ITER"/),&
         description="Maximum number of iterations allowed for the solution of the Z-vector equations.",&
         usage="MAX_ITER  50",&
         default_i_val=30,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="EPS_CONV",&
         description="Convergence threshold for the solution of the Z-vector equations. "//&
                     "The Z-vector equations have the form of a linear system of equations Ax=b, "//&
                     "convergence is achieved when |Ax-b|<=EPS_CONV.",&
         usage="EPS_CONV 1.0E-6",type_of_var=real_t,&
         default_r_val=1.0E-4_dp, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)


  END SUBROUTINE create_cphf

! *****************************************************************************
!> \brief ...
!> \param section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_mp2_potential(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mp2_potential', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,"INTERACTION_POTENTIAL",&
         description="Parameters the interaction potential in computing the biel integrals",&
         n_keywords=1, n_subsections=0, repeats=.FALSE., &
         error=error)

    NULLIFY(keyword)
    CALL keyword_create(&
         keyword=keyword,&
         name="POTENTIAL_TYPE",&
         description="Which interaction potential should be used "//&
                      "(Coulomb, TShPSC operator).",&
         usage="POTENTIAL_TYPE TSHPSC",&
         enum_c_vals=s2a("COULOMB","TShPSC"),&
         enum_i_vals=(/do_hfx_potential_coulomb,&
                       do_mp2_potential_TShPSC/),&
         enum_desc=s2a("Coulomb potential: 1/r", &
                       "TShPSC:<ul><li>1/x - s/Rc for x &le; Rc</li>"//&
                       "<li>(1 - s)/Rc - (x - Rc)/Rc^2 + (x - Rc)^2/Rc^3 - "//&
                       "(2*n^2 - 7*n + 9 - 4*s)*(x - Rc)^3/(Rc^4*(n^2 - 2*n + 1)*(n - 1)) + "//&
                       "(6-3*s - 4*n + n^2)*(x - Rc)^4/(Rc^5*(n^4 - 4*n^3 + 6*n^2 - 4*n + 1)) "//&
                       "for Rc &lt; x &le; n*Rc (4th order polynomial)</li>"//&
                       "<li>0 for x &gt; n*Rc</li></ul>"),&
         default_i_val=do_hfx_potential_coulomb, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="TRUNCATION_RADIUS",&
         description="Determines truncation radius for the truncated TShPSC potential. "//&
                     "Only valid when doing truncated calculation",&
         usage="TRUNCATION_RADIUS 10.0",type_of_var=real_t,&
         default_r_val=10.0_dp,&
         unit_str="angstrom",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(&
         keyword=keyword,&
         name="TShPSC_DATA",&
         description="Location of the file TShPSC.dat that contains the data for the "//&
                     "evaluation of the TShPSC G0 ",&
         usage="TShPSC_DATA t_sh_p_s_c.dat",&
         default_c_val="t_sh_p_s_c.dat",&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)


  END SUBROUTINE create_mp2_potential

END MODULE input_cp2k_mp2

