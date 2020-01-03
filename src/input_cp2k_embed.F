!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief builds the input structure for the EMBED environment: clone of MIXED environment
!> \author Vladimir Rybkin - University of Zurich
! **************************************************************************************************
MODULE input_cp2k_embed
   USE bibliography,                    ONLY: Heaton_Burgess2007,&
                                              Huang2011
   USE cp_output_handling,              ONLY: add_last_numeric,&
                                              cp_print_key_section_create,&
                                              low_print_level
   USE input_constants,                 ONLY: dfet,&
                                              dmfet
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

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_embed'

   PUBLIC :: create_embed_section

CONTAINS

! **************************************************************************************************
!> \brief Create the input section for EMBED: clone of the subroutines for MIXED
!> \param section the section to create
!> \author Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE create_embed_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_embed_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: sub2section, sub3section, subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="EMBED", &
                          description="This section contains all information to run embedded "// &
                          "calculations.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE., &
                          citations=(/Huang2011, Heaton_Burgess2007/))
      NULLIFY (keyword, subsection)

      CALL keyword_create(keyword, __LOCATION__, name="EMBED_METHOD", &
                          description="Select DFET or DMFET.", &
                          usage="EMBED_METHOD DFET", &
                          default_i_val=dfet, &
                          enum_c_vals=s2a("DFET", "DMFET"), &
                          enum_desc=s2a("DFET", "DMFET"), &
                          enum_i_vals=(/dfet, dmfet/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Group partitioning
      CALL keyword_create(keyword, __LOCATION__, name="GROUP_PARTITION", &
                          description="gives the exact number of processors for each group."// &
                          " If not specified processors allocated will be equally distributed for"// &
                          " the specified subforce_eval, trying to build a number of groups equal to the"// &
                          " number of subforce_eval specified.", &
                          usage="group_partition  2 2 4 2 4 ", type_of_var=integer_t, n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NGROUPS", variants=(/"NGROUP"/), &
                          description="Gives the wanted number of groups. Currently must be set to 1", &
                          usage="ngroups 4", type_of_var=integer_t, default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Mapping of atoms
      NULLIFY (sub2section, sub3section)
      CALL section_create(subsection, __LOCATION__, name="MAPPING", &
                          description="Defines the mapping of atoms for the different force_eval with the mixed force_eval."// &
                          " The default is to have a mapping 1-1 between atom index (i.e. all force_eval share the same"// &
                          " geometrical structure). The mapping is based on defining fragments and the mapping the "// &
                          " fragments between the several force_eval and the mixed force_eval", &
                          n_keywords=1, n_subsections=0, repeats=.TRUE.)

      ! Mixed force_eval
      CALL section_create(sub2section, __LOCATION__, name="FORCE_EVAL_EMBED", &
                          description="Defines the fragments for the embedding force_eval (reference)", &
                          n_keywords=1, n_subsections=0, repeats=.TRUE.)

      CALL section_create(sub3section, __LOCATION__, name="FRAGMENT", &
                          description="Fragment definition", &
                          n_keywords=1, n_subsections=0, repeats=.TRUE.)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Defines the index of the fragment defined", &
                          usage="<INTEGER>", type_of_var=integer_t, n_var=1)
      CALL section_add_keyword(sub3section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="_DEFAULT_KEYWORD_", &
                          description="Starting and ending atomic index defining one fragment must be provided", &
                          usage="{Integer} {Integer}", type_of_var=integer_t, n_var=2, repeats=.TRUE.)
      CALL section_add_keyword(sub3section, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(sub2section, sub3section)
      CALL section_release(sub3section)
      CALL section_add_subsection(subsection, sub2section)
      CALL section_release(sub2section)

      ! All other force_eval
      CALL section_create(sub2section, __LOCATION__, name="FORCE_EVAL", &
                          description="Defines the fragments and the mapping for each force_eval (an integer index (ID) "// &
                          "needs to be provided as parameter)", &
                          n_keywords=1, n_subsections=0, repeats=.TRUE.)

      CALL keyword_create( &
         keyword, __LOCATION__, name="DEFINE_FRAGMENTS", &
         description="Specify the fragments definition of the force_eval through the fragments of the"// &
         " force_eval_embed. This avoids the pedantic definition of the fragments for the force_eval,"// &
         " assuming the order of the fragments for the specified force_eval is the same as the sequence "// &
         " of integers provided. Easier to USE should be preferred to the specification of the single fragments.", &
         usage="DEFINE_FRAGMENTS <INTEGER> .. <INTEGER>", type_of_var=integer_t, n_var=-1)
      CALL section_add_keyword(sub2section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Defines the index of the force_eval for which fragments and mappings are provided", &
                          usage="<INTEGER>", type_of_var=integer_t, n_var=1)
      CALL section_add_keyword(sub2section, keyword)
      CALL keyword_release(keyword)

      CALL section_create(sub3section, __LOCATION__, name="FRAGMENT", &
                          description="Fragment definition", &
                          n_keywords=1, n_subsections=0, repeats=.TRUE.)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Defines the index of the fragment defined", &
                          usage="<INTEGER>", type_of_var=integer_t, n_var=1)
      CALL section_add_keyword(sub3section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="_DEFAULT_KEYWORD_", &
                          description="Starting and ending atomic index defining one fragment must be provided", &
                          usage="{Integer} {Integer}", type_of_var=integer_t, n_var=2, repeats=.FALSE.)
      CALL section_add_keyword(sub3section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAP", &
                          description="Provides the index of the fragment of the MIXED force_eval mapped on the"// &
                          " locally defined fragment.", &
                          usage="MAP <INTEGER>", type_of_var=integer_t, n_var=1, repeats=.FALSE.)
      CALL section_add_keyword(sub3section, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(sub2section, sub3section)
      CALL section_release(sub3section)
      CALL section_add_subsection(subsection, sub2section)
      CALL section_release(sub2section)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_print_embed_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)
   END SUBROUTINE create_embed_section

! **************************************************************************************************
!> \brief Create the print section for embedding
!> \param section the section to create
!> \author Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE create_print_embed_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_print_embed_section', &
         routineP = moduleN//':'//routineN

      TYPE(section_type), POINTER                        :: print_key

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="print", &
                          description="Section of possible print options in EMBED env.", &
                          n_keywords=0, n_subsections=1, repeats=.FALSE.)

      NULLIFY (print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "PROGRAM_RUN_INFO", &
                                       description="Controls the printing of information during the evaluation of "// &
                                       "the embedding environment. ", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_print_embed_section

END MODULE input_cp2k_embed
