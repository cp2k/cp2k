!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Functions that build SMEAGOL input section.
!> \author Sergey Chulkov
!> \author Christian Ahart
!> \author Clotilde Cucinotta
! **************************************************************************************************
MODULE input_cp2k_smeagol
   USE bibliography, ONLY: Ahart2024
   USE cp_output_handling, ONLY: add_last_numeric, &
                                 cp_print_key_section_create, &
                                 debug_print_level, &
                                 high_print_level, &
                                 low_print_level, &
                                 medium_print_level, &
                                 silent_print_level
   USE cp_units, ONLY: cp_unit_to_cp2k
   USE input_constants, ONLY: smeagol_bulklead_left, &
                              smeagol_bulklead_leftright, &
                              smeagol_bulklead_right, &
                              smeagol_gridmethod_adaptive, &
                              smeagol_gridmethod_traditional, &
                              smeagol_integraltype_gauss_chebyshev, &
                              smeagol_integraltype_gauss_legendre, &
                              smeagol_runtype_bulktransport, &
                              smeagol_runtype_emtransport
   USE input_keyword_types, ONLY: keyword_create, &
                                  keyword_release, &
                                  keyword_type
   USE input_section_types, ONLY: section_add_keyword, &
                                  section_add_subsection, &
                                  section_create, &
                                  section_release, &
                                  section_type
   USE input_val_types, ONLY: char_t, &
                              integer_t, &
                              logical_t, &
                              real_t
   USE kinds, ONLY: dp
   USE physcon, ONLY: kelvin
   USE string_utilities, ONLY: s2a
#include "./base/base_uses.f90"
   #:include 'input_cp2k_smeagol.fypp'
   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_smeagol'

   PUBLIC :: create_dft_smeagol_section

CONTAINS

! **************************************************************************************************
!> \brief Create SMEAGOL input section.
!> \param section input section
! **************************************************************************************************
   SUBROUTINE create_dft_smeagol_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      ! +11 -- explicitly defined keywords
      !       _SECTION_PARAMETERS_, BS.SubSystemsDelta, AM.AtomListBS, TEMPERATURE, PROJECT_NAME, RUN_TYPE, REGRESSION_TEST,
      !        BulkLead, NIMAGES_IJ, EnergyGridType, TypeOfIntegral
      CALL section_create(section, __LOCATION__, name="SMEAGOL", &
                          description="Parameters for NEGF/SMEAGOL simulation.", &
                          n_subsections=1, repeats=.FALSE., &
                          citations=(/Ahart2024/))

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="_SECTION_PARAMETERS_", &
                          description="Controls activation of SMEAGOL interface", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! +++ reademtr()

      ! BS.SubSystemsBoundaries subsection
      NULLIFY (subsection)
      CALL section_create(subsection, __LOCATION__, name="BS.SubSystemsBoundaries", &
                          description="Bound state subsystems' boundaries", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="_DEFAULT_KEYWORD_", &
                          description="BS.Subsystems lines containing two integer numbers", &
                          usage="lower_bound upper_bound", &
                          repeats=.TRUE., n_var=-1, type_of_var=integer_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! BS.SubSystemsDelta
      CALL keyword_create(keyword, __LOCATION__, name="BS.SubSystemsDelta", &
                          description="Specifies a list of atoms to include.", &
                          usage="BS.SubSystemsDelta {real} {real} .. {real}", repeats=.FALSE., &
                          n_var=-1, type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! use a list of atomic indices (AM.SpeciesBS) instead of list of atomic kind labels (AM.SpeciesBS)
      CALL keyword_create(keyword, __LOCATION__, name="AM.AtomListBS", &
                          description="Specifies a list of atoms to include.", &
                          usage="AM.SpeciesBS {integer} {integer} .. {integer}", repeats=.FALSE., &
                          n_var=-1, type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Electronic temperature of the extended molecule.
      CALL keyword_create(keyword, __LOCATION__, name="TEMPERATURE", &
                          description="Electronic temperature of the extended molecule."// &
                          " Temperature of semi-infinite leads are taken from SMEAR / ELECTRONIC_TEMPERATURE"// &
                          " keyword of relevant bulk-transport calculations.", &
                          n_var=1, type_of_var=real_t, unit_str="K", &
                          default_r_val=300.0_dp/kelvin)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PROJECT_NAME", &
                          description="Name of this SMEAGOL project. SMEAGOL writes certain quantities into files named as "// &
                          "smeagol_project.EXTENSION. The length of the project name should not exceed 20 characters.", &
                          n_var=1, type_of_var=char_t, default_c_val="smeagol_project")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! reademtr() enum keywords

      CALL keyword_create(keyword, __LOCATION__, name="RUN_TYPE", &
                          description="Type of SMEAGOL calculation", &
                          default_i_val=smeagol_runtype_bulktransport, &
                          enum_c_vals=s2a("BulkTransport", "EMTransport"), &
                          enum_desc=s2a( &
                          "Save Kohn-Sham and density matrices of semi-infinite leads.", &
                          "Run transport calculation for the extended molecule."), &
                          enum_i_vals=(/smeagol_runtype_bulktransport, smeagol_runtype_emtransport/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REGRESSION_TEST", &
                          description="Enables regression test mode. "// &
                          "Kohn-Sham and density matrices of semi-infinite leads are not written down.", &
                          n_var=1, type_of_var=logical_t, default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BulkLead", &
                          description="Lead label for bulk transport calculation", &
                          default_i_val=smeagol_bulklead_leftright, &
                          enum_c_vals=s2a("LR", "L", "R"), &
                          enum_desc=s2a( &
                          "Identical left and right leads.", &
                          "Left bulk lead.", &
                          "Right bulk lead."), &
                          enum_i_vals=(/smeagol_bulklead_leftright, smeagol_bulklead_left, smeagol_bulklead_right/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! force the number of cell images along i and j cell vectors
      CALL keyword_create(keyword, __LOCATION__, name="NIMAGES_IJ", &
                          description="Number of cell images along i and j cell vectors. "// &
                          "If missing, determine the number of images based on sparsity pattern of the overlap matrix.", &
                          usage="NIMAGES_IJ 3 3", &
                          n_var=2, type_of_var=integer_t, default_i_vals=(/0, 0/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! reademtr() logical keywords
      #:for name1, keyword1, val1 in reademtr_local_llist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=logical_t, default_l_val=${val1}$, &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      #:for name1, keyword1, val1 in reademtr_negfmod_llist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=logical_t, default_l_val=${val1}$, &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      #:for name1, keyword1, val1 in reademtr_negfcoop_llist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=logical_t, default_l_val=${val1}$, &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      ! reademtr() integer keywords
      #:for name1, keyword1, val1 in reademtr_local_ilist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=integer_t, default_i_val=${val1}$)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      #:for name1, keyword1, val1 in reademtr_negfmod_ilist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=integer_t, default_i_val=${val1}$)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      #:for name1, keyword1, val1 in reademtr_negfcoop_ilist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=integer_t, default_i_val=${val1}$)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      ! reademtr() real-valued keywords
      #:for name1, keyword1, val1 in reademtr_local_rlist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=real_t, default_r_val=${val1}$)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      #:for name1, keyword1, val1 in reademtr_negfmod_rlist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=real_t, default_r_val=${val1}$)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      ! reademtr() physical-valued keywords
      #:for name1, keyword1, val1, unit1 in reademtr_local_plist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="${unit1}$", default_r_val=${val1}$)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      #:for name1, keyword1, val1, unit1 in reademtr_local_explicit_plist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="${unit1}$", default_r_val=${val1}$)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      #:for name1, keyword1, val1, unit1 in reademtr_negfmod_rydberg_plist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="${unit1}$", default_r_val=${val1}$)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      #:for name1, keyword1, val1, unit1 in reademtr_negfmod_plist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="${unit1}$", default_r_val=${val1}$)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      ! +++ ReadOptionsNEGF_DFT()

      ! ReadOptionsNEGF_DFT() logical keywords
      #:for name1, keyword1, val1 in readoptsnegf_negfmod_llist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=logical_t, default_l_val=${val1}$, &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      ! ReadOptionsNEGF_DFT() integer keywords
      #:for name1, keyword1, val1 in readoptsnegf_negfmod_ilist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=integer_t, default_i_val=${val1}$)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      ! ReadOptionsNEGF_DFT() real-valued keywords
      #:for name1, keyword1, val1 in readoptsnegf_negfmod_rlist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=real_t, default_r_val=${val1}$)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      ! ReadOptionsNEGF_DFT() physical-valued keywords
      #:for name1, keyword1, val1, unit1 in readoptsnegf_negfmod_explicit_plist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="${unit1}$", default_r_val=${val1}$)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      #:for name1, keyword1, val1, unit1 in readoptsnegf_bfield_rydberg_plist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="${unit1}$", default_r_val=${val1}$)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      ! +++ emtrans_options()

      ! emtrans_options() enum keywords
      CALL keyword_create(keyword, __LOCATION__, name="EnergyGridType", &
                          description="Type of energy grid", &
                          default_i_val=smeagol_gridmethod_traditional, &
                          enum_c_vals=s2a("TRADITIONAL", "ADAPTIVEGRID"), &
                          enum_desc=s2a( &
                          "Traditional grid.", &
                          "Adaptive grid."), &
                          enum_i_vals=(/smeagol_gridmethod_traditional, smeagol_gridmethod_adaptive/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TypeOfIntegral", &
                          description="Method to integrate Green's functions", &
                          default_i_val=smeagol_integraltype_gauss_legendre, &
                          enum_c_vals=s2a("GAUSS-LEGENDRE", "GAUSS-CHEBYSHEV"), &
                          enum_desc=s2a( &
                          "Gauss-Legendre quadrature method.", &
                          "Gauss-Chebyshev quadrature method."), &
                          enum_i_vals=(/smeagol_integraltype_gauss_legendre, smeagol_integraltype_gauss_chebyshev/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! emtrans_options() logical keywords
      #:for name1, keyword1, val1 in emtoptions_negfmod_llist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=logical_t, default_l_val=${val1}$, &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      #:for name1, keyword1, val1 in emtoptions_local_ilist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=integer_t, default_i_val=${val1}$)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      #:for name1, keyword1, val1 in emtoptions_local_explicit_ilist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=integer_t, default_i_val=${val1}$)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      #:for name1, keyword1, val1 in emtoptions_negfmod_ilist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=integer_t, default_i_val=${val1}$)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      #:for name1, keyword1, val1 in emtoptions_negfmod_explicit_ilist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=integer_t, default_i_val=${val1}$)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      #:for name1, keyword1, val1 in emtoptions_sigma_ilist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=integer_t, default_i_val=${val1}$)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      ! emtrans_options() real-valued keywords
      #:for name1, keyword1, val1 in emtoptions_negfmod_rlist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=real_t, default_r_val=${val1}$)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      #:for name1, keyword1, val1 in emtoptions_negfmod_explicit_rlist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=real_t, default_r_val=${val1}$)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      ! emtrans_options() physical-valued keywords
      #:for name1, keyword1, val1, unit1 in emtoptions_negfmod_rydberg_plist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="${unit1}$", default_r_val=${val1}$)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      ! interface physical-valued keywords
      #:for name1, keyword1, val1, unit1 in smeagol_interface_local_plist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="${unit1}$", default_r_val=${val1}$)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

      #:for name1, keyword1, val1, unit1 in smeagol_interface_local_explicit_plist
         CALL keyword_create(keyword, __LOCATION__, name="${keyword1}$", &
                             description="${keyword1}$", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="${unit1}$", default_r_val=${val1}$)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      #:endfor

   END SUBROUTINE create_dft_smeagol_section

END MODULE input_cp2k_smeagol
