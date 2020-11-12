!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Creates the NNP section of the input
!> \author Christoph Schran (christoph.schran@rub.de)
!> \date   2020-10-10
! **************************************************************************************************
MODULE input_cp2k_nnp

   USE bibliography,                    ONLY: Behler2007,&
                                              Behler2011,&
                                              Schran2020a,&
                                              Schran2020b
   USE cp_output_handling,              ONLY: cp_print_key_section_create,&
                                              medium_print_level
   USE cp_units,                        ONLY: cp_unit_to_cp2k
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_add_subsection,&
                                              section_create,&
                                              section_release,&
                                              section_type
   USE input_val_types,                 ONLY: char_t,&
                                              real_t
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_nnp'

   PUBLIC :: create_nnp_section

CONTAINS

! **************************************************************************************************
!> \brief Create the input section for NNP
!> \param section the section to create
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE create_nnp_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_nnp_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection, subsubsection

!   ------------------------------------------------------------------------

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="NNP", &
                          description="This section contains all information to run a "// &
                          "Neural Network Potential (NNP) calculation.", &
                          n_keywords=3, n_subsections=3, repeats=.FALSE., &
                          citations=(/Behler2007, Behler2011, Schran2020a, Schran2020b/))

      NULLIFY (subsection, subsubsection, keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NNP_INPUT_FILE_NAME", &
                          description="File containing the input information for "// &
                          "the setup of the NNP (n2p2/RuNNer format).", &
                          repeats=.FALSE., default_lc_val="input.nn")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="SCALE_FILE_NAME", &
                          description="File containing the scaling information for "// &
                          "the symmetry functions of the NNP.", &
                          repeats=.FALSE., default_lc_val="scaling.data")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! BIAS subsection
      CALL section_create(subsection, __LOCATION__, name="BIAS", &
                          description="Section to bias the standard deviation of the committee disagreement (sigma) by "// &
                          "E = 0.5 * K_B * (sigma - SIGMA_0)**2, if sigma > SIGMA_0.", &
                          n_keywords=2, n_subsections=0, repeats=.FALSE.)
      CALL keyword_create(keyword, __LOCATION__, name="K_B", &
                          description="Harmonic spring constant of the bias potential [1/hartree].", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.1_dp, &
                          usage="K_B 0.1", &
                          citations=(/Schran2020b/))
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="SIGMA_0", &
                          description="Shift of the harmonic bias potential.", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=cp_unit_to_cp2k(value=0.1_dp, unit_str="hartree"), &
                          unit_str="hartree", &
                          usage="SIGMA_0 [hartree] 0.1")
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="ALIGN_NNP_ENERGIES", &
                          description="Remove PES shifts within the committee by "// &
                          "subtracting energy for each committee member. Provide "// &
                          "one number per C-NNP member.", &
                          repeats=.FALSE., &
                          n_var=-1, &
                          type_of_var=real_t, &
                          usage="ALIGN_NNP_ENERGIES <REAL> <REAL> ... <REAL>")
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      ! print bias subsubsection:
      CALL create_nnp_bias_print_section(subsubsection)
      CALL section_add_subsection(subsection, subsubsection)
      CALL section_release(subsubsection)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)
      ! end BIAS subsection

      CALL section_create(subsection, __LOCATION__, name="MODEL", &
                          description="Section for a single NNP model. "// &
                          "If this section is repeated, a committee model (C-NNP)"// &
                          "is used where the NNP members share the same symmetry functions.", &
                          n_keywords=1, n_subsections=0, repeats=.TRUE.)
      CALL keyword_create(keyword, __LOCATION__, name="WEIGHTS", &
                          description="File containing the weights for the "// &
                          "artificial neural networks of the NNP. "// &
                          "The specified name is extended by .XXX.data", &
                          repeats=.FALSE., default_lc_val="weights")
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_nnp_print_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_nnp_section

! **************************************************************************************************
!> \brief Creates the print section for the nnp subsection
!> \param section the section to create
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE create_nnp_print_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_nnp_print_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key

!   ------------------------------------------------------------------------

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="PRINT", &
                          description="Section of possible print options in NNP code.", &
                          n_keywords=0, n_subsections=5, repeats=.FALSE.)

      NULLIFY (print_key, keyword)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "ENERGIES", &
                                       description="Controls the printing of the NNP energies.", &
                                       print_level=medium_print_level, common_iter_levels=1)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "FORCES", &
                                       description="Controls the printing of the NNP forces.", &
                                       print_level=medium_print_level, common_iter_levels=1)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "FORCES_SIGMA", &
                                       description="Controls the printing of the STD per atom of the NNP forces.", &
                                       print_level=medium_print_level, common_iter_levels=1)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "EXTRAPOLATION", &
                                       description="If activated, output structures with extrapolation "// &
                                       "warning in xyz-format", &
                                       print_level=medium_print_level, common_iter_levels=1)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "SUM_FORCE", &
                                       description="If activated, output summed force over specified atoms. "// &
                                       "Used in Green-Kubo relation for friction at liquid-solid interfaces.", &
                                       print_level=medium_print_level, common_iter_levels=1)

      CALL keyword_create(keyword, __LOCATION__, name="ATOM_LIST", &
                          description="List of atoms over which to calculate summed force", &
                          usage="ATOM_LISTS {O} {H} .. {X}", repeats=.FALSE., &
                          n_var=-1, type_of_var=char_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_nnp_print_section

! **************************************************************************************************
!> \brief Creates the print section for the nnp bias subsubsection
!> \param section the section to create
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE create_nnp_bias_print_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_nnp_bias_print_section', &
         routineP = moduleN//':'//routineN

      TYPE(section_type), POINTER                        :: print_key

!   ------------------------------------------------------------------------

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="PRINT", &
                          description="Section of possible print options in NNP code.", &
                          n_keywords=0, n_subsections=3, repeats=.FALSE.)

      NULLIFY (print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "BIAS_ENERGY", &
                                       description="Controls the printing of the BIAS energy.", &
                                       print_level=medium_print_level, common_iter_levels=1)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "BIAS_FORCES", &
                                       description="Controls the printing of the BIAS forces.", &
                                       print_level=medium_print_level, common_iter_levels=1)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_nnp_bias_print_section

END MODULE input_cp2k_nnp
