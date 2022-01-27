!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief function that build the dft section of the input
!> \par History
!>      10.2005 moved out of input_cp2k [fawzi]
!> \author fawzi
! **************************************************************************************************
MODULE input_cp2k_voronoi
   USE bibliography,                    ONLY: Brehm2018,&
                                              Brehm2020,&
                                              Brehm2021,&
                                              Rycroft2009,&
                                              Thomas2015
   USE cp_output_handling,              ONLY: cp_print_key_section_create,&
                                              debug_print_level
   USE input_constants,                 ONLY: voro_radii_cov,&
                                              voro_radii_unity,&
                                              voro_radii_user,&
                                              voro_radii_vdw
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_type
   USE input_val_types,                 ONLY: integer_t,&
                                              lchar_t,&
                                              real_t
   USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_voronoi'

   PUBLIC :: create_print_voronoi_section

CONTAINS

! **************************************************************************************************
!> \brief Create the print voronoi section
!> \param print_key ...
!> \author Martin Brehm
! **************************************************************************************************
   SUBROUTINE create_print_voronoi_section(print_key)
      TYPE(section_type), POINTER                        :: print_key

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(print_key))

      ! Voronoi Integration via LibVori
      CALL cp_print_key_section_create(print_key, __LOCATION__, name="VORONOI", &
                                       description="Controls the Voronoi integration of the total electron density"// &
                                       " for the computation of electromagnetic moments, see [Thomas2015],"// &
                                       " [Brehm2020], and [Brehm2021] (via LibVori,"// &
                                       " see <a href=""https://brehm-research.de/voronoi"" target=""_blank"">"// &
                                       " https://brehm-research.de/voronoi</a> ).", &
                                       print_level=debug_print_level + 1, filename="", &
                                       citations=(/Rycroft2009, Thomas2015, Brehm2018, Brehm2020, Brehm2021/))

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="Appends frames to already existing .voronoi file.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SANITY_CHECK", &
                          description="Performs a sanity check before each Voronoi integration, i.e.,"// &
                          " checks if every grid point is located in exactly one Voronoi cell.", &
                          usage="SANITY_CHECK T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OVERWRITE", &
                          description="Specify this keyword to overwrite any existing ""properties.emp"" file if"// &
                          " it already exists. By default, the data is appended to an existing .emp file.", &
                          usage="OVERWRITE T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SKIP_FIRST", &
                          description="Skips the first step of a MD run (avoids duplicate step if restarted).", &
                          usage="SKIP_FIRST T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VERBOSE", &
                          description="Switches on verbose screen output of the Voronoi integration.", &
                          usage="VERBOSE T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OUTPUT_EMP", &
                          description="Writes the resulting electromagnetic moments to a binary file ""properties.emp""."// &
                          " The file name cannot be changed.", &
                          usage="OUTPUT_EMP T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OUTPUT_TEXT", &
                          description="Writes the resulting electromagnetic moments to text files (*.voronoi)."// &
                          " The file name is specified via FILENAME.", &
                          usage="OUTPUT_TEXT T", default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REFINEMENT_FACTOR", &
                          description="Sets the refinement factor for the Voronoi integration.", &
                          usage="REFINEMENT 2", n_var=1, default_i_val=1, type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VORONOI_RADII", &
                          description="Which atomic radii to use for the radical Voronoi tessellation.", &
                          usage="VORONOI_RADII {UNITY,VDW,COVALENT,USER}", repeats=.FALSE., n_var=1, &
                          default_i_val=voro_radii_vdw, &
                          enum_c_vals=s2a("UNITY", "VDW", "COVALENT", "USER"), &
                          enum_desc=s2a("Use unity radii (non-radical Voronoi tessellation)", "Use VdW atom radii", &
                                        "Use covalent atom radii", "Use user-specified atom radii"), &
                          enum_i_vals=(/voro_radii_unity, voro_radii_vdw, voro_radii_cov, voro_radii_user/))
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="USER_RADII", &
                          description="Defines user atom radii for the radical Voronoi tessellation (one per atom).", &
                          usage="USER_RADII {real} {real} {real}", repeats=.FALSE., &
                          unit_str="angstrom", &
                          type_of_var=real_t, n_var=-1)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MOLECULAR_PROPERTIES", &
                          description="Calculation of molecular properties from Voronoi integration.", &
                          usage="MOLECULAR_PROPERTIES T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MOLPROP_FILE_NAME", &
                          description="Root of the file name where to print molecular properties."// &
                          " filename.molprop is used.", &
                          usage="MOLPROP_FILE_NAME <FILENAME>", &
                          default_lc_val="__STD_OUT__", type_of_var=lchar_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_print_voronoi_section

END MODULE input_cp2k_voronoi
