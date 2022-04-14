!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief input section for MP2
!> \par History
!>      05.2011 created
!> \author MDB
! **************************************************************************************************
MODULE input_cp2k_mp2
   USE bibliography,                    ONLY: &
        Bates2013, DelBen2012, DelBen2013, DelBen2015, DelBen2015b, Rybkin2016, Wilhelm2016a, &
        Wilhelm2016b, Wilhelm2017, Wilhelm2018
   USE cp_eri_mme_interface,            ONLY: create_eri_mme_section
   USE cp_output_handling,              ONLY: add_last_numeric,&
                                              cp_print_key_section_create,&
                                              debug_print_level,&
                                              high_print_level,&
                                              low_print_level,&
                                              medium_print_level,&
                                              silent_print_level
   USE cp_units,                        ONLY: cp_unit_to_cp2k
   USE input_constants,                 ONLY: &
        do_eri_gpw, do_eri_mme, do_eri_os, do_potential_coulomb, do_potential_id, &
        do_potential_long, do_potential_mix_cl, do_potential_short, do_potential_truncated, do_potential_tshpsc, &
        eri_default, gaussian, gw_no_print_exx, gw_pade_approx, gw_print_exx, gw_read_exx, &
        gw_skip_for_regtest, gw_two_pole_model, kp_weights_W_auto, kp_weights_W_tailored, &
        kp_weights_W_uniform, mp2_method_direct, mp2_method_gpw, mp2_method_none, numerical, &
        ri_default, ri_rpa_g0w0_crossing_bisection, ri_rpa_g0w0_crossing_newton, &
        ri_rpa_g0w0_crossing_z_shot, wfc_mm_style_gemm, wfc_mm_style_syrk, &
        gw_gf_mic, gw_gf_gamma, ot_precond_full_all, ot_precond_full_kinetic, ot_precond_full_single, &
      ot_precond_full_single_inverse, ot_precond_none, ot_precond_s_inverse
   USE input_cp2k_hfx,                  ONLY: create_hfx_section
   USE input_cp2k_kpoints,              ONLY: create_kpoint_set_section
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
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_mp2'

   PUBLIC :: create_mp2_section

CONTAINS

! **************************************************************************************************
!> \brief creates the input section for the mp2 part
!> \param section the section to create
!> \author MDB
! **************************************************************************************************
   SUBROUTINE create_mp2_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="WF_CORRELATION", &
                          description="Sets up the wavefunction-based correlation methods as MP2, "// &
                          "RI-MP2, RI-SOS-MP2, RI-RPA and GW (inside RI-RPA). ", &
                          n_keywords=4, n_subsections=7, repeats=.TRUE., &
                          citations=(/DelBen2012, DelBen2013, DelBen2015, DelBen2015b, Rybkin2016, &
                                      Wilhelm2016a, Wilhelm2016b, Wilhelm2017, Wilhelm2018/))

      NULLIFY (keyword, subsection)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="MEMORY", &
         description="Maximum allowed total memory usage during MP2 methods [MiB].", &
         usage="MEMORY 1500 ", &
         default_r_val=1.024E+3_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="SCALE_S", &
         description="Scaling factor of the singlet energy component (opposite spin, OS) of the "// &
         "MP2, RI-MP2 and SOS-MP2 correlation energy. ", &
         usage="SCALE_S  1.0", &
         default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="SCALE_T", &
         description="Scaling factor of the triplet energy component (same spin, SS) of the MP2 "// &
         "and RI-MP2 correlation energy.", &
         usage="SCALE_T  1.0", &
         default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="GROUP_SIZE", &
         variants=(/"NUMBER_PROC"/), &
         description="Group size used in the computation of GPW and MME integrals and the MP2 correlation energy. "// &
         "A smaller GROUP_SIZE (for example GROUP_SIZE = number of MPI ranks per node) "// &
         "accelerates the computation of integrals but a too small "// &
         "GROUP_SIZE may lead to out of memory. To use all processors set GROUP_SIZE -1 "// &
         "(might lead to VERY high computation times). GROUP_SIZE must be a divisor of the total number "// &
         "of MPI ranks.", &
         usage="GROUP_SIZE 32", &
         default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (subsection)
      CALL create_mp2_details_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_ri_mp2(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_ri_rpa(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_ri_laplace(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! here we generate an imag. time subsection to use with RPA or Laplace-SOS-MP2
      CALL create_low_scaling(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_ri_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_integrals_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (print_key)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "PRINT", &
                                       description="Controls the printing basic info about WFC methods", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_mp2_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_mp2_details_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="MP2", &
                          description="Parameters influencing MP2 (non-RI).", &
                          n_keywords=3, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Activates MP2 calculations.", &
                          usage="&MP2 .TRUE.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="METHOD", &
         citations=(/DelBen2012, DelBen2013/), &
         description="Method that is used to compute the MP2 energy.", &
         usage="METHOD MP2_GPW", &
         enum_c_vals=s2a("NONE", "DIRECT_CANONICAL", "MP2_GPW"), &
         enum_i_vals=(/mp2_method_none, mp2_method_direct, mp2_method_gpw/), &
         enum_desc=s2a("Skip MP2 calculation.", &
                       "Use the direct mp2 canonical approach.", &
                       "Use the GPW approach to MP2."), &
         default_i_val=mp2_method_direct)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="BIG_SEND", &
         description="Influencing the direct canonical MP2 method: Send big "// &
         "messages between processes (useful for >48 processors).", &
         usage="BIG_SEND", &
         default_l_val=.TRUE., &
         lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_mp2_details_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_ri_mp2(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="RI_MP2", &
                          description="Parameters influencing the RI-MP2 method. RI-MP2 supports gradients.", &
                          n_keywords=3, n_subsections=1, repeats=.FALSE., &
                          citations=(/DelBen2013/))

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Putting the &RI_MP2 section activates RI-MP2 calculation.", &
                          usage="&RI_MP2 .TRUE.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (subsection)
      CALL create_cphf(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL keyword_create(keyword, __LOCATION__, name="BLOCK_SIZE", &
                          variants=(/"MESSAGE_SIZE"/), &
                          description="Determines the blocking used for communication in RI-MP2. Larger BLOCK_SIZE "// &
                          "reduces communication but requires more memory. The default (-1) is automatic.", &
                          usage="BLOCK_SIZE 2", &
                          default_i_val=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_CANONICAL", &
                          description="Threshold for discriminate if a given ij pairs "// &
                          "of the unrelaxed MP2 density matrix has to be "// &
                          "calculated with a canonical reformulation based "// &
                          "on the occupied eigenvalues differences.", &
                          usage="EPS_CANONICAL 1.0E-8", type_of_var=real_t, &
                          default_r_val=1.0E-7_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="FREE_HFX_BUFFER", &
         description="Free the buffer containing the 4 center integrals used in the Hartree-Fock exchange calculation. "// &
         "This will be effective only for gradients calculations, since for the energy only "// &
         "case, the buffers are released by default. (Right now debugging only).", &
         usage="FREE_HFX_BUFFER", &
         default_l_val=.TRUE., &
         lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="USE_OLD_GRADIENT_CODE", &
         description="Use the original RI-MP2 gradient code.", &
         usage="USE_OLD_GRADIENT_CODE  T", &
         lone_keyword_l_val=.TRUE., &
         default_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="PRINT_DGEMM_INFO", &
         description="Print details about all DGEMM calls.", &
         lone_keyword_l_val=.TRUE., &
         default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_ri_mp2

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_opt_ri_basis(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="OPT_RI_BASIS", &
                          description="Parameters influencing the optimization of the RI MP2 basis. "// &
                          "Only exponents of non-contracted auxiliary basis can be optimized. "// &
                          "An initial RI auxiliary basis has to be specified.", &
                          n_keywords=6, n_subsections=0, repeats=.FALSE., &
                          citations=(/DelBen2013/))
      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Putting the &OPT_RI_BASIS section activates optimization of RI basis.", &
                          usage="&OPT_RI_BASIS .TRUE.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DELTA_I_REL", &
                          variants=(/"DI_REL"/), &
                          description="Target accuracy in the relative deviation of the amplitudes calculated with "// &
                          "and without RI approximation, (more details in Chem.Phys.Lett.294(1998)143).", &
                          usage="DELTA_I_REL  1.0E-6_dp", &
                          default_r_val=1.0E-6_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DELTA_RI", &
                          variants=(/"DRI"/), &
                          description="Target accuracy in the absolute difference between the RI-MP2 "// &
                          "and the exact MP2 energy, DRI=ABS(E_MP2-E_RI-MP2).", &
                          usage="DELTA_RI  1.0E-6_dp", &
                          default_r_val=5.0E-6_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_DERIV", &
                          variants=(/"EPS_NUM_DERIV"/), &
                          description="The derivatives of the MP2 energy with respect to the "// &
                          "exponents of the basis are calculated numerically. "// &
                          "The change in the exponent a_i employed for the numerical evaluation "// &
                          "is defined as h_i=EPS_DERIV*a_i.", &
                          usage="EPS_DERIV  1.0E-3_dp", &
                          default_r_val=1.0E-3_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER", &
                          variants=(/"MAX_NUM_ITER"/), &
                          description="Specifies the maximum number of steps in the RI basis optimization.", &
                          usage="MAX_ITER 100", &
                          default_i_val=50)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NUM_FUNC", &
                          description="Specifies the number of function, for each angular momentum (s, p, d ...), "// &
                          "employed in the automatically generated initial guess. "// &
                          "This will be effective only if RI_AUX_BASIS_SET in the KIND section is not specified.", &
                          usage="NUM_FUNC {number of s func.} {number of p func.} ...", &
                          n_var=-1, default_i_vals=(/-1/), type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BASIS_SIZE", &
                          description="Specifies the size of the auxiliary basis set automatically "// &
                          "generated as initial guess. This will be effective only if RI_AUX_BASIS_SET "// &
                          "in the KIND section and NUM_FUNC are not specified.", &
                          usage="BASIS_SIZE  (MEDIUM|LARGE|VERY_LARGE)", &
                          enum_c_vals=s2a("MEDIUM", "LARGE", "VERY_LARGE"), &
                          enum_i_vals=(/0, 1, 2/), &
                          default_i_val=0)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_opt_ri_basis

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_ri_laplace(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="RI_SOS_MP2", &
                          description="Parameters influencing the RI-SOS-MP2-Laplace method", &
                          n_keywords=3, n_subsections=1, repeats=.FALSE., &
                          citations=(/DelBen2013/))

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Putting the &RI_SOS_MP2 section activates RI-SOS-MP2 calculation.", &
                          usage="&RI_SOS_MP2 .TRUE.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="QUADRATURE_POINTS", &
         variants=(/"LAPLACE_NUM_QUAD_POINTS"/), &
         description="Number of quadrature points for the numerical integration in the RI-SOS-MP2-Laplace method.", &
         usage="QUADRATURE_POINTS 6", &
         default_i_val=5)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="SIZE_INTEG_GROUP", &
         variants=(/"LAPLACE_GROUP_SIZE"/), &
         description="Group size for the integration in the Laplace method, that is the number of processes involved in "// &
         "the computation of each integration point. SIZE_INTEG_GROUP has to be a multiple "// &
         "of GROUP_SIZE in the WF_CORRELATION section. The default (-1) "// &
         "is automatic.", &
         usage="SIZE_INTEG_GROUP 16", &
         default_i_val=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_ri_laplace

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_ri_rpa(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="RI_RPA", &
                          description="Parameters influencing RI-RPA and GW.", &
                          n_keywords=8, n_subsections=4, repeats=.FALSE., &
                          citations=(/DelBen2013, DelBen2015/))

      NULLIFY (keyword, subsection)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Putting the &RI_RPA section activates RI-RPA calculation.", &
                          usage="&RI_RPA .TRUE.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="QUADRATURE_POINTS", &
                          variants=(/"RPA_NUM_QUAD_POINTS"/), &
                          description="Number of quadrature points for the numerical integration in the RI-RPA method.", &
                          usage="QUADRATURE_POINTS 60", &
                          default_i_val=40)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SIZE_FREQ_INTEG_GROUP", &
                          variants=(/"RPA_GROUP_SIZE"/), &
                          description="Group size for frequency integration, that is the number of processes involved in "// &
                          "the computation of each integration point. SIZE_FREQ_INTEG_GROUP has to be a multiple "// &
                          "of GROUP_SIZE in the WF_CORRELATION section. The default (-1) "// &
                          "is automatic.", &
                          usage="SIZE_FREQ_INTEG_GROUP 16", &
                          default_i_val=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="MM_STYLE", &
                          description="Matrix multiplication style for the Q matrix.", &
                          usage="MM_STYLE GEMM", &
                          enum_c_vals=s2a("GEMM", "SYRK"), &
                          enum_i_vals=(/wfc_mm_style_gemm, wfc_mm_style_syrk/), &
                          enum_desc=s2a("Use pdgemm: more flops, maybe faster.", &
                                        "Use pdysrk: fewer flops, maybe slower."), &
                          default_i_val=wfc_mm_style_gemm)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="MINIMAX_QUADRATURE", &
         variants=(/"MINIMAX"/), &
         description="Use the Minimax quadrature scheme for performing the numerical integration. "// &
         "Maximum number of quadrature point limited to 20.", &
         usage="MINIMAX_QUADRATURE", &
         default_l_val=.FALSE., &
         lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="RI_AXK", &
         variants=(/"AXK"/), &
         description="Decide whether to perform an RPA-AXK calculation.", &
         usage="RI_AXK", &
         default_l_val=.FALSE., &
         lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="RSE", &
         variants=(/"SE"/), &
         description="Decide whether to add singles correction.", &
         usage="RI_AXK", &
         default_l_val=.FALSE., &
         lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="ADMM", &
         description="Decide whether to perform ADMM in the exact exchange calc. for RPA and/or GW. "// &
         "Not really recommended. In most cases, the Hartree-Fock exchange is not expensive "// &
         "and there is no need for ADMM. Only in case of diffuse functions for periodic systems "// &
         "(only reasonable for GW bandgap calculations) exact exchange can be expensive and then, "// &
         "we recommend using RI for exchange.", &
         usage="ADMM", &
         default_l_val=.FALSE., &
         lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="SCALE_RPA", &
         description="Scales RPA energy contributions (RPA, AXK).", &
         usage="SCALE_RPA 1.0", &
         default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="PRINT_DGEMM_INFO", &
         description="Print details about all DGEMM calls.", &
         lone_keyword_l_val=.TRUE., &
         default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! here we generate a hfx subsection to use in the case EXX has to be computed after RPA
      CALL create_hfx_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! here we generate a G0W0 subsection to use if G0W0 is desired
      CALL create_ri_g0w0(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! here we generate a RI_AXK subsection
      CALL create_ri_axk(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_ri_rpa

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_ri_axk(section)
      TYPE(section_type), POINTER                        :: section

      !TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="RI_AXK", &
                          description="Parameters influencing the RI-RPA-AXK method", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE., &
                          citations=(/Bates2013/))

      !NULLIFY (keyword, subsection)

   END SUBROUTINE create_ri_axk

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_ri_g0w0(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="GW", &
                          description="Parameters influencing the RI-G0W0 method", &
                          n_keywords=24, n_subsections=1, repeats=.FALSE.)

      NULLIFY (keyword, subsection)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Activates cubic-scaling RPA, GW and Laplace-SOS-MP2 calculations.", &
                          usage="&GW .TRUE.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CORR_MOS_OCC", &
                          variants=(/"CORR_OCC"/), &
                          description="Number of occupied MOs whose energies are corrected by RI-G0W0. "// &
                          "Counting beginning from HOMO, e.g. 3 corrected occ. MOs correspond "// &
                          "to correction of HOMO, HOMO-1 and HOMO-2. Numerical effort and "// &
                          "storage of RI-G0W0 increase linearly with this number. In case you "// &
                          "want to correct all occ. MOs, insert a number larger than the number "// &
                          "of occ. MOs.", &
                          usage="CORR_OCC 3", &
                          default_i_val=10)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CORR_MOS_VIRT", &
                          variants=(/"CORR_VIRT"/), &
                          description="Number of virtual MOs whose energies are corrected by RI-G0W0. "// &
                          "Counting beginning from LUMO, e.g. 3 corrected occ. MOs correspond "// &
                          "to correction of LUMO, LUMO+1 and LUMO+2. Numerical effort and "// &
                          "storage of RI-G0W0 increase linearly with this number. In case you "// &
                          "want to correct all virt. MOs, insert a number larger than the number "// &
                          "of virt. MOs.", &
                          usage="CORR_VIRT 3", &
                          default_i_val=10)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NUMB_POLES", &
                          description="Number of poles for the fitting. Usually, two poles are sufficient. ", &
                          usage="NUMB_POLES 2", &
                          default_i_val=2)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OMEGA_MAX_FIT", &
                          description="Determines fitting range for the self-energy on the imaginary axis: "// &
                          "[0, OMEGA_MAX_FIT] for virt orbitals, [-OMEGA_MAX_FIT,0] for occ orbitals. "// &
                          "Unit: Hartree. Default: 0.734996 H = 20 eV. ", &
                          usage="OMEGA_MAX_FIT 0.5", &
                          default_r_val=0.734996_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CROSSING_SEARCH", &
                          description="Determines, how the self_energy is evaluated on the real axis.", &
                          usage="CROSSING_SEARCH Z_SHOT", &
                          enum_c_vals=s2a("Z_SHOT", "NEWTON", "BISECTION"), &
                          enum_i_vals=(/ri_rpa_g0w0_crossing_z_shot, &
                                        ri_rpa_g0w0_crossing_newton, ri_rpa_g0w0_crossing_bisection/), &
                          enum_desc=s2a("Calculate the derivative of Sigma and out of it Z. Then extrapolate using Z.", &
                                        "Make a Newton-Raphson fix point iteration.", &
                                        "Make a bisection fix point iteration."), &
                          default_i_val=ri_rpa_g0w0_crossing_newton)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FERMI_LEVEL_OFFSET", &
                          description="Fermi level for occ. orbitals: e_HOMO + FERMI_LEVEL_OFFSET; "// &
                          "Fermi level for virt. orbitals: e_LUMO - FERMI_LEVEL_OFFSET. "// &
                          "In case e_homo + FERMI_LEVEL_OFFSET < e_lumo - FERMI_LEVEL_OFFSET, "// &
                          "we set Fermi level = (e_HOMO+e_LUMO)/2. For cubic-scaling GW, the Fermi level "// &
                          "is always equal to (e_HOMO+e_LUMO)/2 regardless of FERMI_LEVEL_OFFSET.", &
                          usage="FERMI_LEVEL_OFFSET 1.0E-2", &
                          default_r_val=2.0E-2_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EV_GW_ITER", &
                          description="Maximum number of iterations for eigenvalue self-consistency cycle. The "// &
                          "computational effort of GW scales linearly with this number. ", &
                          usage="EV_GW_ITER 3", &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SC_GW0_ITER", &
                          description="Maximum number of iterations for GW0 self-consistency cycle. The "// &
                          "computational effort of GW is not much affected by the number of scGW0 cycles. ", &
                          usage="SC_GW0_ITER 3", &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_ITER", &
                          description="Target accuracy for the eigenvalue self-consistency. "// &
                          "If the G0W0 HOMO-LUMO gap differs by less than the "// &
                          "target accuracy during the iteration, the eigenvalue "// &
                          "self-consistency cycle stops. Unit: Hartree.", &
                          usage="EPS_EV_SC_ITER 0.00005", &
                          default_r_val=cp_unit_to_cp2k(value=0.00136_dp, unit_str="eV"), &
                          unit_str="eV")

      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PRINT_EXX", &
                          description="Print exchange self-energy minus exchange correlation potential for Gamma-only "// &
                          "calculation (PRINT). For a GW calculation with k-points we use this output as "// &
                          "exchange self-energy (READ). This is a temporary solution because the hybrid MPI/OMP "// &
                          "parallelization in the HFX by Manuel Guidon conflicts with the parallelization in "// &
                          "low-scaling GW k-points which is most efficient with maximum number of MPI tasks and "// &
                          "minimum number of OMP threads. For HFX by M. Guidon, the density matrix is "// &
                          "fully replicated on every MPI rank which necessitates a high number of OMP threads per MPI "// &
                          "rank for large systems to prevent out of memory. "// &
                          "Such a high number of OMP threads would slow down the GW calculation "// &
                          "severely. Therefore, it was decided to temporarily divide the GW k-point calculation in a "// &
                          "Gamma-only HF calculation with high number of OMP threads to prevent out of memory and "// &
                          "a GW k-point calculation with 1 OMP thread per MPI rank reading the previousHF output.", &
                          usage="PRINT_EXX TRUE", &
                          enum_c_vals=s2a("TRUE", "FALSE", "READ", "SKIP_FOR_REGTEST"), &
                          enum_i_vals=(/gw_print_exx, gw_no_print_exx, gw_read_exx, gw_skip_for_regtest/), &
                          enum_desc=s2a("Please, put TRUE for Gamma only calculation to get the exchange self-energy. "// &
                                        "If 'SIGMA_X' and the corresponding values for the exchange-energy are written, "// &
                                        "the writing has been successful", &
                                        "FALSE is needed if you want to do nothing here.", &
                                        "Please, put READ for the k-point GW calculation to read the exact exchange. "// &
                                        "You have to provide an output file including the exact exchange. This file "// &
                                        "has to be named 'exx.dat'.", &
                                        "SKIP_FOR_REGTEST is only used for the GW k-point regtest where no exchange "// &
                                        "self-energy is computed."), &
                          default_i_val=gw_no_print_exx)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RI_SIGMA_X", &
                          description="If true, the exchange self-energy is calculated approximatively with RI. "// &
                          "If false, the Hartree-Fock implementation in CP2K is used.", &
                          usage="RI_SIGMA_X", &
                          default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="IC_CORR_LIST", &
                          description="List of image charge correction from a previous calculation to be applied in G0W0 "// &
                          "or evGW. Keyword is active, if the first entry is positive (since IC corrections are positive "// &
                          "occupied MOs. The start corresponds to the first corrected GW level.", &
                          usage="IC_CORR_LIST <REAL> ... <REAL>", &
                          default_r_vals=(/-1.0_dp/), &
                          type_of_var=real_t, n_var=-1, unit_str="eV")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="IC_CORR_LIST_BETA", &
                          description="IC_CORR_LIST for beta spins in case of open shell calculation.", &
                          usage="IC_CORR_LIST_BETA <REAL> ... <REAL>", &
                          default_r_vals=(/-1.0_dp/), &
                          type_of_var=real_t, n_var=-1, unit_str="eV")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PERIODIC_CORRECTION", &
                          description="If true, the periodic correction scheme is used employing k-points. "// &
                          "Method is not recommended to use, use instead PERIODIC_LOW_SCALING which much "// &
                          "more accurate than the periodic correction.", &
                          usage="PERIODIC_CORRECTION", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BSE", &
                          description="If true, electronic excitation energies are computed from the Bethe-"// &
                          "Salpeter equation on top of GW eigenvalues.  Parameter of BSE can be adjusted in "// &
                          "the corresponding section.", &
                          usage="BSE", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="IMAGE_CHARGE_MODEL", &
                          variants=(/"IC"/), &
                          description="If true, an image charge model is applied to mimic the renormalization of "// &
                          "electronic levels of a molecule at a metallic surface. For this calculation, the molecule "// &
                          "has to be reflected on the desired xy image plane. The coordinates of the reflected molecule "// &
                          "have to be added to the coord file as ghost atoms. For the ghost atoms, identical basis sets "// &
                          "the normal atoms have to be used.", &
                          usage="IC TRUE", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ANALYTIC_CONTINUATION", &
                          description="Defines which type of analytic continuation for the self energy is used", &
                          usage="ANALYTIC_CONTINUATION", &
                          enum_c_vals=s2a("TWO_POLE", "PADE"), &
                          enum_i_vals=(/gw_two_pole_model, gw_pade_approx/), &
                          enum_desc=s2a("Use 'two-pole' model.", &
                                        "Use Pade approximation."), &
                          default_i_val=gw_pade_approx)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NPARAM_PADE", &
                          description="Number of parameters for the Pade approximation "// &
                          "when using the latter for the analytic continuation of the "// &
                          "self energy. 16 parameters (corresponding to 8 poles) are "// &
                          "are recommended.", &
                          usage="NPARAM_PADE 16", &
                          default_i_val=16)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GAMMA_ONLY_SIGMA", &
                          variants=(/"GAMMA"/), &
                          description="If true, the correlation self-energy is only computed at the Gamma point. "// &
                          "The Gamma point itself is obtained by averaging over all kpoints of the DFT mesh.", &
                          usage="GAMMA TRUE", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="UPDATE_XC_ENERGY", &
                          description="If true, the Hartree-Fock and RPA total energy are printed and the total energy "// &
                          "is corrected using exact exchange and the RPA correlation energy.", &
                          usage="UPDATE_XC_ENERGY", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="KPOINTS_SELF_ENERGY", &
                          description="Specify number of k-points for the k-point grid of the self-energy. Internally, a "// &
                          "Monkhorst-Pack grid is used. A dense k-point grid may be necessary to compute an accurate density "// &
                          "of state from GW. Large self-energy k-meshes do not cost much more computation time.", &
                          usage="KPOINTS  nx  ny  nz", repeats=.TRUE., &
                          n_var=3, type_of_var=integer_t, default_i_vals=(/0, 0, 0/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! here we generate a subsection for the periodic GW correction
      CALL create_periodic_gw_correction_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! here we generate a subsection for Bethe-Salpeter
      CALL create_bse_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! here we generate a subsection for image charge calculations
      CALL create_ic_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! here we generate a subsection for additional printing
      CALL create_print_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! here we generate a subsection for calculating the GW band structures
      CALL create_kpoint_set_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_ri_g0w0

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_print_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="PRINT", &
                          description="Parameters for GW related output.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="", &
                          variants=(/"SELF_ENERGY"/), &
                          description="If true, print the self-energy for all levels for real energy "// &
                          "together with the straight line to see the quasiparticle energy as intersection. "// &
                          "In addition, prints the self-energy for imaginary frequencies together with the Pade fit.", &
                          usage="SELF_ENERGY TRUE", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_periodic_gw_correction_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="PERIODIC_CORRECTION", &
                          description="Parameters influencing correction for periodic GW. Old method, "// &
                          "not recommended to use", &
                          n_keywords=12, n_subsections=1, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="KPOINTS", &
                          description="Specify number of k-points for a single k-point grid. Internally, a "// &
                          "Monkhorst-Pack grid is used. Typically, even numbers are chosen such that the Gamma "// &
                          "point is excluded from the k-point mesh.", &
                          usage="KPOINTS  nx  ny  nz", repeats=.TRUE., &
                          n_var=3, type_of_var=integer_t, default_i_vals=(/16, 16, 16/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NUM_KP_GRIDS", &
                          description="Number of k-point grids around the Gamma point with different resolution. "// &
                          "E.g. for KPOINTS 4 4 4 and NUM_KP_GRIDS 3, there will be a 3x3x3 Monkhorst-Pack (MP) k-point "// &
                          "grid for the whole Brillouin zone (excluding Gamma), another 3x3x3 MP grid with smaller  "// &
                          "spacing around Gamma (again excluding Gamma) and a very fine 4x4x4 MP grid around Gamma.", &
                          usage="NUM_KP_GRIDS 5", &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_KPOINT", &
                          description="If the absolute value of a k-point is below EPS_KPOINT, this kpoint is "// &
                          "neglected since the Gamma point is not included in the periodic correction.", &
                          usage="EPS_KPOINT 1.0E-4", &
                          default_r_val=1.0E-05_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MO_COEFF_GAMMA", &
                          description="If true, only the MO coefficients at the Gamma point are used for the periodic  "// &
                          "correction. Otherwise, the MO coeffs are computed at every k-point which is much more "// &
                          "expensive. It should be okay to use the Gamma MO coefficients.", &
                          usage="MO_COEFF_GAMMA", &
                          default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="AVERAGE_DEGENERATE_LEVELS", &
                          variants=(/"ADL"/), &
                          description="If true, the correlation self-energy of degenerate levels is averaged.", &
                          usage="AVERAGE_DEGENERATE_LEVELS", &
                          default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_EIGENVAL", &
                          description="Threshold for considering levels as degenerate. Unit: Hartree.", &
                          usage="EPS_EIGENVAL 1.0E-5", &
                          default_r_val=2.0E-04_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EXTRAPOLATE_KPOINTS", &
                          variants=(/"EXTRAPOLATE"/), &
                          description="If true, extrapolates the k-point mesh. Only working if k-point mesh numbers are "// &
                          "divisible by 4, e.g. 8x8x8 or 12x12x12 is recommended.", &
                          usage="EXTRAPOLATE_KPOINTS FALSE", &
                          default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DO_AUX_BAS_GW", &
                          description="If true, use a different basis for the periodic correction. This can be necessary "// &
                          "in case a diffused basis is used for GW to converge the HOMO-LUMO gap. In this case, "// &
                          "numerical problems may occur due to diffuse functions in the basis. This keyword only works if "// &
                          "AUX_GW <basis set>  is specified in the kind section for every atom kind.", &
                          usage="AUX_BAS_GW TRUE", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FRACTION_AUX_MOS", &
                          description="Fraction how many MOs are used in the auxiliary basis.", &
                          usage="FRACTION_AUX_MOS 0.6", &
                          default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NUM_OMEGA_POINTS", &
                          description="Number of Clenshaw-Curtis integration points for the periodic correction in cubic- "// &
                          "scaling GW. This variable is a dummy variable for canonical N^4 GW calculations.", &
                          usage="NUM_OMEGA_POINTS 200", &
                          default_i_val=300)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_bse_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="BSE", &
                          description="Parameters influencing the Bethe-Salpeter calculations on "// &
                          "top of GW eigenvalues.", &
                          n_keywords=12, n_subsections=1, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NUM_Z_VECTORS", &
                          description="Number Z vectors used in the subspace iterations. This is a convergence "// &
                          "parameter: Increasing NUM_Z_VECTORS should keep the result constant in case of convergence,.", &
                          usage="NUM_Z_VECTORS 50", &
                          default_i_val=20)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="THRESHOLD_MIN_TRANS", &
                          variants=(/"EPS"/), &
                          description="Threshold for stopping the iteration for computing the transition energies. "// &
                          "If the lowest excitation changes by less than THRESHOLD_MIN_TRANS (in eV), the ieration"// &
                          "stops.", &
                          usage="THRESHOLD_MIN_TRANS 0.001", &
                          default_r_val=0.001_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER", &
                          description="Maximum number of iterations for determining the transition energies.", &
                          usage="MAX_ITER 200", &
                          default_i_val=200)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_ic_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="IC", &
                          description="Parameters influencing the image charge correction. "// &
                          "The image plane is always an xy plane, so adjust the molecule according "// &
                          "to that. ", &
                          n_keywords=3, n_subsections=1, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PRINT_IC_LIST", &
                          description="If true, the image charge correction values are printed in a list, "// &
                          "such that it can be used as input for a subsequent evGW calculation.", &
                          usage="PRINT_IC_VALUES", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_DIST", &
                          description="Threshold where molecule and image molecule have to coincide. ", &
                          usage="EPS_DIST 0.1", unit_str="angstrom", &
                          type_of_var=real_t, default_r_val=3.0E-02_dp, repeats=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_low_scaling(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create( &
         section, __LOCATION__, name="LOW_SCALING", &
         description="Cubic scaling RI-RPA, GW and Laplace-SOS-MP2 method using the imaginary time formalism. "// &
         "EPS_GRID in WFC_GPW section controls accuracy / req. memory for 3-center integrals. "// &
         "SORT_BASIS EXP should be specified in DFT section.", &
         n_keywords=12, n_subsections=2, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Activates cubic-scaling RPA, GW and Laplace-SOS-MP2 calculations.", &
                          usage="&LOW_SCALING .TRUE.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MEMORY_CUT", &
                          description="Reduces memory for sparse tensor contractions by this factor. "// &
                          "A high value leads to some loss of performance. "// &
                          "This memory reduction factor applies to storage of the tensors 'M occ' / 'M virt' "// &
                          "but does not reduce storage of '3c ints'.", &
                          usage="MEMORY_CUT 16", &
                          default_i_val=5)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MEMORY_INFO", &
                          description="Decide whether to print memory info on the sparse matrices.", &
                          usage="MEMORY_INFO", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="EPS_FILTER", &
         description="Determines a threshold for the DBCSR based multiply."// &
         "Normally, this EPS_FILTER determines accuracy and timing of low-scaling RPA and GW calculations.", &
         usage="EPS_FILTER 1.0E-10 ", type_of_var=real_t, &
         default_r_val=1.0E-9_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="EPS_FILTER_FACTOR", &
         description="Multiply EPS_FILTER with this factor to determine filter epsilon "// &
         "for DBCSR based multiply P(it)=(Mocc(it))^T*Mvirt(it) "// &
         "Default should be kept.", &
         type_of_var=real_t, &
         default_r_val=10.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="EPS_STORAGE_SCALING", &
         variants=(/"EPS_STORAGE"/), &
         description="Scaling factor to scale EPS_FILTER. Storage threshold for compression "// &
         "will be EPS_FILTER*EPS_STORAGE_SCALING.", &
         default_r_val=1.0E-3_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="DO_KPOINTS", &
         description="Besides in DFT, this keyword has to be switched on if one wants to do kpoints in. "// &
         "cubic RPA.", &
         usage="DO_KPOINTS", &
         default_l_val=.FALSE., &
         lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="KPOINTS", &
         description="Keyword activates periodic, low-scaling GW calculations (&LOW_SCALING section also needed). "// &
         "For periodic calculations, kpoints are used for the density response, the "// &
         "Coulomb interaction and the screened Coulomb interaction. For 2d periodic systems, e.g. xz "// &
         "periodicity, please also specify KPOINTS, e.g.  N_x  1  N_z.", &
         usage="KPOINTS  N_x  N_y  N_z", &
         n_var=3, type_of_var=integer_t, default_i_vals=(/0, 0, 0/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="KPOINT_WEIGHTS_W", &
         description="For kpoints in low-scaling GW, a Monkhorst-Pack mesh is used. The screened Coulomb "// &
         "interaction W(k) needs special care near the Gamma point (e.g. in 3d, W(k) diverges at the "// &
         "Gamma point with W(k) ~ k^alpha). KPOINT_WEIGHTS_W decides how the weights of the "// &
         "Monkhorst-Pack mesh are chosen to compute W(R) = int_BZ W(k) exp(ikR) dk (BZ=Brllouin zone). ", &
         usage="KPOINT_WEIGHTS_W AUTO", &
         enum_c_vals=s2a("TAILORED", "AUTO", "UNIFORM"), &
         enum_i_vals=(/kp_weights_W_tailored, kp_weights_W_auto, kp_weights_W_uniform/), &
         enum_desc=s2a("Choose k-point integration weights such that the function f(k)=k^alpha is "// &
                       "exactly integrated. alpha is specified using EXPONENT_TAILORED_WEIGHTS.", &
                       "As 'TAILORED', but alpha is chosen automatically according to dimensionality "// &
                       "(3D: alpha = -2 for 3D, 2D: alpha = -1 for exchange self-energy, uniform "// &
                       "weights for correlation self-energy).", &
                       "Choose the same weight for every k-point (original Monkhorst-Pack method)."), &
         default_i_val=kp_weights_W_auto)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="EXPONENT_TAILORED_WEIGHTS", &
         description="Gives the exponent of exactly integrated function in case 'KPOINT_WEIGHTS_W "// &
         "TAILORED' is chosen.", &
         usage="EXPONENT_TAILORED_WEIGHTS -2", &
         default_r_val=-2.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="REL_CUTOFFS_CHI_W", &
         description="For periodic calculations, the irreducible density reponse chi^0(r,r') and the screened "// &
         "Coulomb interaction W(r,r') is truncated "// &
         "in real space. No truncation of chi^0(r,r') and W(r,r') for |r-r'| < cutoff_1*(smallest lattice vector length), "// &
         "chi^0(r,r') is set to zero for |r-r'| > cutoff_2*(smallest lattice vector length). Smooth transition "// &
         "in between. For cutoff_1 = 0.5, and cutoff_2 = 0.5, we have the minimum image convention (that is also default).", &
         usage="REL_CUTOFFS_CHI cutoff_1 cutoff_2", &
         n_var=2, type_of_var=real_t, default_r_vals=(/0.5_dp, 0.5_dp/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="REGULARIZATION_RI", &
         description="Parameter to reduce the expansion coefficients in RI for periodic GW. Larger parameter "// &
         "means smaller expansion coefficients that leads to a more stable calculation at the price "// &
         "of a slightly worse RI approximation. In case the parameter 0.0 is chosen, ordinary RI is used.", &
         usage="REGULARIZATION_RI 1.0E-4", &
         default_r_val=1.0E-2_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="MAKE_CHI_POS_DEFINITE", &
         description="If true, makes eigenvalue decomposition of chi(iw,k) and removes negative "// &
         "eigenvalues. May increase computational cost significantly. Only recommended to try in case "// &
         "Cholesky decomposition of epsilon(iw,k) fails.", &
         usage="MAKE_CHI_POS_DEFINITE", &
         default_l_val=.TRUE., &
         lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="TRUNC_COLOUMB_RI_X", &
         description="If true, use the truncated Coulomb operator for the exchange-self-energy in "// &
         "periodic GW.",&
         usage="TRUNCATED_COLOUMB_RI_X", &
         default_l_val=.FALSE., &
         lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="KEEP_QUADRATURE", &
         variants=s2a("KEEP_WEIGHTS", "KEEP_QUAD", "KEEP_WEIGHT"), &
         description="Keep the Laplace quadrature defined at the first energy evaluations throughout "// &
         "the run. Allows to have consistent force evaluations.", &
         usage="KEEP_QUADRATURE", &
         default_l_val=.TRUE., &
         lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="K_MESH_G_FACTOR", &
         description="The k-mesh for the Green's function can be chosen to be larger than the k-mesh for "// &
         "W (without much higher computational cost). The factor given here multiplies the mesh for W to obtain"// &
         "the k-mesh for G. Example: factor 4, k-mesh for W: 4x4x1 -> k-mesh for G: 16x16x1 (z-dir. is "// &
         "non-periodic).", &
         default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (keyword)
      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="GREENS_FUNCTION", &
         description="Decides which scheme is used for the periodic Green's function "// &
                     "(only relevant for periodic GW calculations)."// &
                     "Both schemes are exact in the limit of large unit cells."  , &
         usage="GREENS_FUNCTION GAMMA", &
         enum_c_vals=s2a("GAMMA", "MIC"), &
         enum_i_vals=(/gw_gf_gamma, gw_gf_mic/), &
         enum_desc=s2a("Use the Green's functions at the Gamma point. ", &
                       "Use the minimum-image convention for the Green's function."), &
         default_i_val=gw_gf_gamma)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="MIN_BLOCK_SIZE", &
         description="Minimum tensor block size. Adjusting this value may have minor effect on "// &
         "performance but default should be good enough.", &
         default_i_val=5)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="MIN_BLOCK_SIZE_MO", &
         description="Tensor block size for MOs. Only relevant for GW calculations. "// &
         "The memory consumption of GW scales as O(MIN_BLOCK_SIZE_MO). It is recommended to "// &
         "set this parameter to a smaller number if GW runs out of memory. "// &
         "Otherwise the default should not be changed.", &
         default_i_val=64)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (subsection)
      CALL create_low_scaling_cphf(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_wfc_gpw(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="WFC_GPW", &
                          description="Parameters for the GPW approach in Wavefunction-based Correlation methods", &
                          n_keywords=5, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="EPS_GRID", &
                          description="Determines a threshold for the GPW based integration", &
                          usage="EPS_GRID 1.0E-9 ", type_of_var=real_t, &
                          default_r_val=1.0E-8_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="EPS_FILTER", &
         description="Determines a threshold for the DBCSR based multiply (usually 10 times smaller than EPS_GRID)."// &
         "Normally, this EPS_FILTER determines accuracy and timing of cubic-scaling RPA calculation.", &
         usage="EPS_FILTER 1.0E-10 ", type_of_var=real_t, &
         default_r_val=1.0E-9_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF", &
                          description="The cutoff of the finest grid level in the MP2 gpw integration.", &
                          usage="CUTOFF 300", type_of_var=real_t, &
                          default_r_val=300.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REL_CUTOFF", &
                          variants=(/"RELATIVE_CUTOFF"/), &
                          description="Determines the grid at which a Gaussian is mapped.", &
                          usage="REL_CUTOFF 50", type_of_var=real_t, &
                          default_r_val=50.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PRINT_LEVEL", &
                          variants=(/"IOLEVEL"/), &
                          description="How much output is written by the individual groups.", &
                          usage="PRINT_LEVEL HIGH", &
                          default_i_val=silent_print_level, enum_c_vals= &
                          s2a("SILENT", "LOW", "MEDIUM", "HIGH", "DEBUG"), &
                          enum_desc=s2a("Almost no output", &
                                        "Little output", "Quite some output", "Lots of output", &
                                        "Everything is written out, useful for debugging purposes only"), &
                          enum_i_vals=(/silent_print_level, low_print_level, medium_print_level, &
                                        high_print_level, debug_print_level/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="EPS_PGF_ORB_S", &
         description="Screening for overlap matrix in RI. Usually, it is best to choose this parameter "// &
         "to be very small since the inversion of overlap matrix might be ill-conditioned.", &
         usage="EPS_PGF_ORB_S 1.0E-10 ", type_of_var=real_t, &
         default_r_val=1.0E-10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_wfc_gpw

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_cphf(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create( &
         section, __LOCATION__, name="CPHF", &
         description="Parameters influencing the solution of the Z-vector equations in MP2 gradients calculations.", &
         n_keywords=2, n_subsections=0, repeats=.FALSE., &
         citations=(/DelBen2013/))

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER", &
                          variants=(/"MAX_NUM_ITER"/), &
                          description="Maximum number of iterations allowed for the solution of the Z-vector equations.", &
                          usage="MAX_ITER  50", &
                          default_i_val=30)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_CONV", &
                          description="Convergence threshold for the solution of the Z-vector equations. "// &
                          "The Z-vector equations have the form of a linear system of equations Ax=b, "// &
                          "convergence is achieved when |Ax-b|<=EPS_CONV.", &
                          usage="EPS_CONV 1.0E-6", type_of_var=real_t, &
                          default_r_val=1.0E-4_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_cphf

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_low_scaling_cphf(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="CPHF", &
                          description="Parameters influencing the solution of the Z-vector equations "// &
                          "in low-scaling Laplace-SOS-MP2 gradients calculations.", &
                          n_keywords=5, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_CONV", &
                          description="Target accuracy for Z-vector euation solution.", &
                          usage="EPS 1.e-6", default_r_val=1.e-6_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER", &
                          description="Maximum number of conjugate gradient iteration to be performed for one optimization.", &
                          usage="MAX_ITER 200", default_i_val=50)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="PRECONDITIONER", &
         description="Type of preconditioner to be used with all minimization schemes. "// &
         "They differ in effectiveness, cost of construction, cost of application. "// &
         "Properly preconditioned minimization can be orders of magnitude faster than doing nothing.", &
         usage="PRECONDITIONER FULL_ALL", &
         default_i_val=ot_precond_full_all, &
         enum_c_vals=s2a("FULL_ALL", "FULL_SINGLE_INVERSE", "FULL_SINGLE", "FULL_KINETIC", "FULL_S_INVERSE", &
                         "NONE"), &
         enum_desc=s2a("Most effective state selective preconditioner based on diagonalization, "// &
                       "requires the ENERGY_GAP parameter to be an underestimate of the HOMO-LUMO gap. "// &
                       "This preconditioner is recommended for almost all systems, except very large systems where "// &
                       "make_preconditioner would dominate the total computational cost.", &
                       "Based on H-eS cholesky inversion, similar to FULL_SINGLE in preconditioning efficiency "// &
                       "but cheaper to construct, "// &
                       "might be somewhat less robust. Recommended for large systems.", &
                       "Based on H-eS diagonalisation, not as good as FULL_ALL, but somewhat cheaper to apply. ", &
                       "Cholesky inversion of S and T, fast construction, robust, and relatively good, "// &
                       "use for very large systems.", &
                       "Cholesky inversion of S, not as good as FULL_KINETIC, yet equally expensive.", &
                       "skip preconditioning"), &
         enum_i_vals=(/ot_precond_full_all, ot_precond_full_single_inverse, ot_precond_full_single, &
                       ot_precond_full_kinetic, ot_precond_s_inverse, ot_precond_none/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ENERGY_GAP", &
                          description="Energy gap estimate [a.u.] for preconditioning", &
                          usage="ENERGY_GAP 0.1", &
                          default_r_val=0.2_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_low_scaling_cphf

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_mp2_potential(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="INTERACTION_POTENTIAL", &
                          description="Parameters the interaction potential in computing the biel integrals", &
                          n_keywords=4, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="POTENTIAL_TYPE", &
         description="Which interaction potential should be used "// &
         "(Coulomb, TShPSC operator).", &
         usage="POTENTIAL_TYPE TSHPSC", &
         enum_c_vals=s2a("COULOMB", "TShPSC", "LONGRANGE", "SHORTRANGE", "TRUNCATED", "MIX_CL"), &
         enum_i_vals=(/do_potential_coulomb, &
                       do_potential_TShPSC, &
                       do_potential_long, &
                       do_potential_short, &
                       do_potential_truncated, &
                       do_potential_mix_cl/), &
         enum_desc=s2a("Coulomb potential: 1/r", &
                       "TShPSC:<ul><li>1/x - s/Rc for x &le; Rc</li>"// &
                       "<li>(1 - s)/Rc - (x - Rc)/Rc^2 + (x - Rc)^2/Rc^3 - "// &
                       "(2*n^2 - 7*n + 9 - 4*s)*(x - Rc)^3/(Rc^4*(n^2 - 2*n + 1)*(n - 1)) + "// &
                       "(6-3*s - 4*n + n^2)*(x - Rc)^4/(Rc^5*(n^4 - 4*n^3 + 6*n^2 - 4*n + 1)) "// &
                       "for Rc &lt; x &le; n*Rc (4th order polynomial)</li>"// &
                       "<li>0 for x &gt; n*Rc</li></ul>", &
                       "Longrange Coulomb potential: erf(w*r)/r)", &
                       "Shortrange Coulomb potential: erfc(w*r)/r", &
                       "Truncated Coulomb potential", &
                       "Mixed Coulomb/Longrange Coulomb potential"), &
         default_i_val=do_potential_coulomb)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TRUNCATION_RADIUS", &
                          variants=(/"CUTOFF_RADIUS"/), &
                          description="Determines truncation radius for the truncated potentials. "// &
                          "Only valid when doing truncated calculations", &
                          usage="TRUNCATION_RADIUS 10.0", type_of_var=real_t, &
                          default_r_val=10.0_dp, &
                          unit_str="angstrom")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="POTENTIAL_DATA", &
         variants=s2a("TShPSC_DATA", "T_C_G_DATA"), &
         description="Location of the file TShPSC.dat or t_c_g.dat that contains the data for the "// &
         "evaluation of the evaluation of the truncated potentials", &
         usage="TShPSC_DATA t_sh_p_s_c.dat", &
         default_c_val="t_sh_p_s_c.dat")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="OMEGA", &
         description="Range separation parameter for the longrange or shortrange potential. "// &
         "Only valid when longrange or shortrange potential is requested.", &
         usage="OMEGA 0.5", type_of_var=real_t, &
         default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="SCALE_COULOMB", &
         description="Scaling factor of (truncated) Coulomb potential in mixed (truncated) Coulomb/Longrange potential."// &
         "Only valid when mixed potential is requested.", &
         usage="OMEGA 0.5", type_of_var=real_t, &
         default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="SCALE_LONGRANGE", &
         description="Scaling factor of longrange Coulomb potential in mixed (truncated) Coulomb/Longrange potential."// &
         "Only valid when mixed potential is requested.", &
         usage="OMEGA 0.5", type_of_var=real_t, &
         default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_mp2_potential

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_ri_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="RI", &
                          description="Parameters influencing resolution of the identity (RI) that is"// &
                          "used in RI-MP2, RI-RPA, RI-SOS-MP2 and GW (inside RI-RPA).", &
                          n_keywords=6, n_subsections=2, repeats=.FALSE.)

      NULLIFY (subsection)
      CALL create_RI_metric_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_opt_ri_basis(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (keyword)
      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="ROW_BLOCK", &
         variants=(/"ROW_BLOCK_SIZE"/), &
         description="Size of the row block used in the SCALAPACK block cyclic data distribution."// &
         "Default is (ROW_BLOCK=-1) is automatic. "// &
         "A proper choice can speedup the parallel matrix multiplication in the case of RI-RPA and RI-SOS-MP2-Laplace.", &
         usage="ROW_BLOCK 512", &
         default_i_val=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="COL_BLOCK", &
         variants=(/"COL_BLOCK_SIZE"/), &
         description="Size of the column block used in the SCALAPACK block cyclic data distribution."// &
         "Default is (COL_BLOCK=-1) is automatic. "// &
         "A proper choice can speedup the parallel matrix multiplication in the case of RI-RPA and RI-SOS-MP2-Laplace.", &
         usage="COL_BLOCK 512", &
         default_i_val=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="CALC_COND_NUM", &
         variants=(/"CALC_CONDITION_NUMBER"/), &
         description="Calculate the condition number of the (P|Q) matrix for the RI methods.", &
         usage="CALC_COND_NUM", &
         default_l_val=.FALSE., &
         lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DO_SVD", &
                          description="Wether to perform a singular value decomposition instead of the Cholesky decomposition "// &
                          "of the potential operator in the RI basis. Computationally expensive but numerically more stable. "// &
                          "It reduces the computational costs of some subsequent steps. Recommended when a longrange Coulomb "// &
                          "potential is employed.", &
                          usage="DO_SVD  .TRUE.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ERI_BLKSIZE", &
                          description="block sizes for tensors (only used if ERI_METHOD=MME). First value "// &
                          "is the block size for ORB basis, second value is the block size for RI_AUX basis.", &
                          usage="ERI_BLKSIZE", &
                          n_var=2, &
                          default_i_vals=(/4, 16/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_ri_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_integrals_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="INTEGRALS", &
                          description="Parameters controlling how to compute integrals that are needed "// &
                          "in MP2, RI-MP2, RI-RPA, RI-SOS-MP2 and GW (inside RI-RPA).", &
                          n_keywords=2, n_subsections=3, repeats=.FALSE.)

      NULLIFY (subsection)
      CALL create_eri_mme_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_wfc_gpw(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_mp2_potential(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="ERI_METHOD", &
                          description="Method for calculating periodic electron repulsion integrals "// &
                          "(MME method is faster but experimental, forces not yet implemented)."// &
                          "Obara-Saika (OS) for the Coulomb operator can only be used for non-periodic calculations.", &
                          usage="ERI_METHOD MME", &
                          enum_c_vals=s2a("DEFAULT", "GPW", "MME", "OS"), &
                          enum_i_vals=(/eri_default, do_eri_gpw, do_eri_mme, do_eri_os/), &
                          enum_desc=s2a("Use default ERI method (for periodic systems: GPW, for molecules: OS, "// &
                                        "for MP2 and RI-MP2: GPW in any case).", &
                                        "Uses Gaussian Plane Wave method [DelBen2013].", &
                                        "Uses MiniMax-Ewald method (experimental, ERI_MME subsection, only for fully periodic "// &
                                        "systems with orthorhombic cells).", &
                                        "Use analytical Obara-Saika method."), &
                          default_i_val=eri_default)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SIZE_LATTICE_SUM", &
                          description="Size of sum range L. ", &
                          usage="SIZE_LATTICE_SUM  10", &
                          default_i_val=5)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_integrals_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_RI_metric_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="RI_METRIC", &
                          description="Sets up RI metric", &
                          repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="POTENTIAL_TYPE", &
         description="Decides which operator/metric is used for resolution of the identity (RI).", &
         usage="POTENTIAL_TYPE DEFAULT", &
         enum_c_vals=s2a("DEFAULT", "COULOMB", "IDENTITY", "LONGRANGE", "SHORTRANGE", "TRUNCATED"), &
         enum_i_vals=(/ri_default, do_potential_coulomb, do_potential_id, do_potential_long, &
                       do_potential_short, do_potential_truncated/), &
         enum_desc=s2a("Use Coulomb metric for RI-MP2 and normal-scaling RI-SOS-MP2, RI-RPA and GW. "// &
                       "Use Overlap metric for low-scaling RI-SOS-MP2, RI-RPA and GW for periodic systems. "// &
                       "Use truncated Coulomb metric for low-scaling RI-SOS-MP2, RI-RPA and GW for non-periodic systems.", &
                       "Coulomb metric: 1/r. Recommended for RI-MP2,", &
                       "Overlap metric: delta(r).", &
                       "Longrange metric: erf(omega*r)/r. Not recommended with DO_SVD .TRUE.", &
                       "Shortrange metric: erfc(omega*r)/r", &
                       "Truncated Coulomb metric: if (r &lt; R_c) 1/r else 0. More "// &
                       "accurate than IDENTITY for non-periodic systems. Recommended for low-scaling methods."), &
         default_i_val=ri_default)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (keyword)
      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="OMEGA", &
         description="The range parameter for the short/long range operator (in 1/a0).", &
         usage="OMEGA 0.5", &
         default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF_RADIUS", &
                          description="The cutoff radius (in Angstrom) for the truncated Coulomb operator.", &
                          usage="CUTOFF_RADIUS 3.0", default_r_val=cp_unit_to_cp2k(value=3.0_dp, unit_str="angstrom"), &
                          type_of_var=real_t, unit_str="angstrom")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="T_C_G_DATA", &
         description="Location of the file t_c_g.dat that contains the data for the "// &
         "evaluation of the truncated gamma function ", &
         default_c_val="t_c_g.dat")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_RANGE", &
                          description="The threshold to determine the effective range of the short range "// &
                          "RI metric: erfc(omega*eff_range)/eff_range = EPS_RANGE", &
                          default_r_val=1.0E-08_dp, &
                          repeats=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_RI_metric_section

END MODULE input_cp2k_mp2
