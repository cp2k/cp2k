!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief function that build the projection of MO in RTP section of the input
!> \author Guillaume Le Breton 04.2023
! **************************************************************************************************
MODULE input_cp2k_projection_rtp
   USE cp_output_handling,              ONLY: cp_print_key_section_create
   USE input_constants,                 ONLY: proj_mo_ref_scf,&
                                              proj_mo_ref_xas_tdp
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_add_subsection,&
                                              section_create,&
                                              section_release,&
                                              section_type
   USE input_val_types,                 ONLY: integer_t
   USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_projection_rtp'

   PUBLIC :: create_projection_rtp_section

CONTAINS

! **************************************************************************************************
!> \brief creates the section for time dependent projection of the MOs
!> \param section ...
!> \author Guillaume Le Breton
! **************************************************************************************************
   SUBROUTINE create_projection_rtp_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="PROJECTION_MO", &
                          description="Projects the Time Dependent (TD) MO "// &
                          "coefficients to reference ones. You can define "// &
                          "several sections like this to project the TD-MOs "// &
                          "on different reference MOs. Note that each section "// &
                          "projects from one spin of the TD MOs (TD_MO_INDEX) to "// &
                          "one spin of the reference MOs (REF_MO_INDEX).", &
                          n_keywords=7, n_subsections=1, &
                          repeats=.TRUE.)

      NULLIFY (keyword, subsection)
      CALL keyword_create(keyword, __LOCATION__, name="REFERENCE_TYPE", &
                          description="Type of the reference MO file provided in REF_MO_FILE_NAME.", &
                          enum_c_vals=s2a("SCF", "XAS_TDP"), &
                          usage="REFERENCE_TYPE SCF", &
                          default_i_val=proj_mo_ref_scf, &
                          enum_desc=s2a("The reference MO is from an SCF calculation.", &
                                        "The reference MO is from an XAS_TDP analysis."), &
                          enum_i_vals=(/proj_mo_ref_scf, proj_mo_ref_xas_tdp/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REF_MO_FILE_NAME", &
                          description="Name of the wavefunction file to read the reference MO from. "// &
                          "For instance, a restart wfn file from SCF calculation or an excited state from XAS_TDP calculation. "// &
                          "If no file is specified, the default is to use DFT%WFN_RESTART_FILE_NAME. "// &
                          "Currently, a RTP restart file (.rtpwfn) cannot be used as reference. "// &
                          "Currently, this file SHALL have the same number of spin as the propagated one "// &
                          "(eventhough you use only the first spin from this reference).", &
                          usage="REF_MO_FILE_NAME <FILENAME>", &
                          default_lc_val="DEFAULT")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REF_MO_INDEX", &
                          description="Indexes of the reference MO read from the .wfn reference file (see REF_MO_FILE_NAME). "// &
                          "Use this keywork if REFERENCE_TYPE=SCF. "// &
                          "Set to -1 to project on all the MO available. "// &
                          "One file will be generated per index defined.", &
                          usage="MO_REF_INDEX 1 2", &
                          default_i_vals=(/1/), &
                          n_var=-1, type_of_var=integer_t, repeats=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REF_MO_SPIN", &
                          description="Spin of the reference MOs to consider. "// &
                          "1 for ALPHA and 2 for BETA spin respectively. "// &
                          "If the reference MO is spin independent this key is not used.", &
                          usage="REF_MO_SPIN 1", &
                          default_i_val=1, &
                          n_var=1, type_of_var=integer_t, repeats=.FALSE.)

      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SUM_ON_ALL_REF", &
            description="Set to .TRUE. in order to sum all the projection done over the required MO_REF_INDEX for each TD MOs. "// &
                          "Only one file will be generated containing the results for every MO_TD_INDEX.", &
                          usage="SUM_ON_ALL_REF .TRUE.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE., &
                          repeats=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TD_MO_INDEX", &
                          description="Indexes of the time dependent MOs to project on the reference MOs. "// &
                          "Set to -1 to project on all the TD MOs.", &
                          usage="MO_TD_INDEX 1 2", &
                          default_i_vals=(/1/), &
                          n_var=-1, type_of_var=integer_t, repeats=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TD_MO_SPIN", &
                          description="Spin of the TD MOs to consider. 1 for ALPHA spin, 2 for BETA spin. "// &
                          "If the TD calculation is spin independent this key is not used.", &
                          usage="MO_TD_SPIN 1", &
                          default_i_val=1, &
                          n_var=1, type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SUM_ON_ALL_TD", &
               description="Set to .TRUE. in order to sum the projection done over all on TD MOs on the required MO_REF_INDEX. "// &
                          "One file per MO_REF_INDEX will be generated. "// &
                          "Combining SUM_ON_ALL_TD and SUM_ON_ALL_REF lead to one file one projection: "// &
                          "the population of all the defined TD_MO_INDEX over the reference MO_REF_INDEX per time step required.", &
                          usage="SUM_ON_ALL_TD .TRUE.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE., &
                          repeats=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! The results for different time step are stored in the same file by default:
      CALL cp_print_key_section_create(subsection, __LOCATION__, name="PRINT", &
                                       description="How to print the MO projection", &
                                       common_iter_levels=999999999, &
                                       filename="PROJ_MO")
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_projection_rtp_section

END MODULE input_cp2k_projection_rtp
