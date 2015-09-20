!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief builds the input structure for the MIXED environment
!> \par History
!>      10.2008 created [tlaino]
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
MODULE input_cp2k_mixed
  USE cp_output_handling,              ONLY: add_last_numeric,&
                                             cp_print_key_section_create,&
                                             low_print_level,&
                                             medium_print_level
  USE input_constants,                 ONLY: mix_coupled,&
                                             mix_generic,&
                                             mix_linear_combination,&
                                             mix_minimum,&
                                             mix_restrained
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: char_t,&
                                             integer_t,&
                                             lchar_t,&
                                             real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_mixed'

  PUBLIC :: create_mix_section

CONTAINS

! *****************************************************************************
!> \brief Create the input section for MIXED.
!> \param section the section to create
!> \author fschiff
! *****************************************************************************
  SUBROUTINE create_mix_section(section)
    TYPE(section_type), POINTER              :: section

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mix_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: sub2section, sub3section, &
                                                subsection

    failure=.FALSE.

    CPASSERT(.NOT.ASSOCIATED(section))
    CALL section_create(section,name="MIXED",&
         description="This section contains all information to run with a hamiltonian "//&
                     "defined by a mixing of force_evals",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)
    NULLIFY(keyword, subsection)

    CALL keyword_create(keyword, name="MIXING_TYPE",&
         description="The type of mixing to be employed",&
         usage="MIXING_TYPE LINEAR_COMBINATION",&
         default_i_val=mix_linear_combination,&
         enum_c_vals=s2a("LINEAR_COMBINATION", &
                         "MINIMUM",&
                         "COUPLED",&
                         "RESTRAINT",&
                         "GENMIX"),&
         enum_desc=s2a("Linear combination of force envs (support only 2 force_evals)", &
                       "Use the force env with the minimum energy (support only 2 force_evals)",&
                       "Consider the force envs as a two state system with a given"//&
                       " coupling matrix element (support only 2 force_evals)",&
                       "Use the difference between the energy of the force envs as a"//&
                       " restraint on the first (support only 2 force_evals)",&
                       "Defines a user-driven generica coupling (support for an unlimited number of force_eval)"),&
         enum_i_vals=(/mix_linear_combination,mix_minimum,mix_coupled,mix_restrained,mix_generic/))
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="GROUP_PARTITION",&
         description="gives the exact number of processors for each group."//&
         " If not specified processors allocated will be equally distributed for"//&
         " the specified subforce_eval, trying to build a number of groups equal to the"//&
         " number of subforce_eval specified.",&
         usage="group_partition  2 2 4 2 4 ", type_of_var=integer_t, n_var=-1)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="NGROUPS",variants=(/"NGROUP"/),&
         description="Gives the wanted number of groups. If not specified the number"//&
         " of groups is set to the number of subforce_eval defined.",&
         usage="ngroups 4", type_of_var=integer_t)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    ! Double force_eval
    CALL section_create(subsection,name="LINEAR",&
         description="Linear combination between two force_eval:  F= lambda F1 + (1-lambda) F2",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)
    CALL keyword_create(keyword, name="LAMBDA",&
         description="Specify the mixing parameter lambda in the formula.",&
         usage="lambda <REAL>", type_of_var=real_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)
    CALL section_add_subsection(section,subsection)
    CALL section_release(subsection)

    CALL section_create(subsection,name="COUPLING",&
         description="Coupling between two force_eval: E=(E1+E2 - sqrt((E1-E2)**2+4*H12**2))/2",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)
    CALL keyword_create(keyword, name="COUPLING_PARAMETER",&
         description="Coupling parameter H12 used in the coupling",&
         usage="COUPLING_PARAMETER <REAL>", type_of_var=real_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)
    CALL section_add_subsection(section,subsection)
    CALL section_release(subsection)

    CALL section_create(subsection,name="RESTRAINT",&
         description="Restraint between two force_eval: E = E1 + k*(E1-E2-t)**2",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)
    CALL keyword_create(keyword, name="RESTRAINT_TARGET",&
         description="Target value of the restraint (t) ",&
         usage="RESTRAINT_TARGET <REAL>", type_of_var=real_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="RESTRAINT_STRENGTH",&
         description="Strength of the restraint (k) in "//&
         "k*(E1-E2-t)**2" ,&
         usage="RESTRAINT_STRENGTH <REAL>", type_of_var=real_t)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)
    CALL section_add_subsection(section,subsection)
    CALL section_release(subsection)

    ! Multiple force_eval
    CALL section_create(subsection,name="GENERIC",&
         description="User driven coupling between two or more force_eval.",&
         n_keywords=1, n_subsections=0, repeats=.FALSE.)
    CALL keyword_create(keyword, name="MIXING_FUNCTION",&
         description="Specifies the mixing functional form in mathematical notation.",&
         usage="MIXING_FUNCTION (E1+E2-LOG(E1/E2))",  type_of_var=lchar_t,&
         n_var=1)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="VARIABLES",&
         description="Defines the variables of the functional form. To allow an efficient"//&
         " mapping the order of the energy variables will be considered identical to the"//&
         " order of the force_eval in the force_eval_order list.",&
         usage="VARIABLES x",  type_of_var=char_t,&
         n_var=-1)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="PARAMETERS",&
         description="Defines the parameters of the functional form",&
         usage="PARAMETERS a b D",  type_of_var=char_t,&
         n_var=-1, repeats=.TRUE.)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="VALUES",&
         description="Defines the values of parameter of the functional form",&
         usage="VALUES ",  type_of_var=real_t,&
         n_var=-1, repeats=.TRUE., unit_str="internal_cp2k")
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="UNITS",&
         description="Optionally, allows to define valid CP2K unit strings for each parameter value. "//&
                     "It is assumed that the corresponding parameter value is specified in this unit.",&
         usage="UNITS angstrom eV*angstrom^-1 angstrom^1 K",  type_of_var=char_t,&
         n_var=-1, repeats=.TRUE.)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="DX",&
         description="Parameter used for computing the derivative with the Ridders' method.",&
         usage="DX <REAL>", default_r_val=0.1_dp, unit_str="bohr")
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="ERROR_LIMIT",&
         description="Checks that the error in computing the derivative is not larger than "//&
         "the value set; in case error is larger a warning message is printed.",&
         usage="ERROR_LIMIT <REAL>", default_r_val=1.0E-12_dp)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)
    CALL section_add_subsection(section,subsection)
    CALL section_release(subsection)

    ! Mapping of atoms
    NULLIFY(sub2section, sub3section)
    CALL section_create(subsection,name="MAPPING",&
         description="Defines the mapping of atoms for the different force_eval with the mixed force_eval."//&
         " The default is to have a mapping 1-1 between atom index (i.e. all force_eval share the same"//&
         " geometrical structure). The mapping is based on defining fragments and the mapping the "//&
         " fragments between the several force_eval and the mixed force_eval",&
         n_keywords=1, n_subsections=0, repeats=.TRUE.)

    ! Mixed force_eval
    CALL section_create(sub2section,name="FORCE_EVAL_MIXED",&
         description="Defines the fragments for the mixed force_eval (reference)",&
         n_keywords=1, n_subsections=0, repeats=.TRUE.)

    CALL section_create(sub3section,name="FRAGMENT",&
         description="Fragment definition",&
         n_keywords=1, n_subsections=0, repeats=.TRUE.)

    CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
         description="Defines the index of the fragment defined",&
         usage="<INTEGER>",  type_of_var=integer_t, n_var=1)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
         description="Starting and ending atomic index defining one fragment must be provided",&
         usage="{Integer} {Integer}", type_of_var=integer_t, n_var=2, repeats=.TRUE.)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

    CALL section_add_subsection(sub2section,sub3section)
    CALL section_release(sub3section)
    CALL section_add_subsection(subsection,sub2section)
    CALL section_release(sub2section)

    ! All other force_eval
    CALL section_create(sub2section,name="FORCE_EVAL",&
         description="Defines the fragments and the mapping for each force_eval (an integer index (ID) "//&
         "needs to be provided as parameter)",&
         n_keywords=1, n_subsections=0, repeats=.TRUE.)

    CALL keyword_create(keyword, name="DEFINE_FRAGMENTS",&
         description="Specify the fragments definition of the force_eval through the fragments of the"//&
         " force_eval_mixed. This avoids the pedantic definition of the fragments for the force_eval,"//&
         " assuming the order of the fragments for the specified force_eval is the same as the sequence "//&
         " of integers provided. Easier to USE should be preferred to the specification of the single fragments.",&
         usage="DEFINE_FRAGMENTS <INTEGER> .. <INTEGER>", type_of_var=integer_t, n_var=-1)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
         description="Defines the index of the force_eval for which fragments and mappings are provided",&
         usage="<INTEGER>",  type_of_var=integer_t, n_var=1)
    CALL section_add_keyword(sub2section,keyword)
    CALL keyword_release(keyword)

    CALL section_create(sub3section,name="FRAGMENT",&
         description="Fragment definition",&
         n_keywords=1, n_subsections=0, repeats=.TRUE.)

    CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
         description="Defines the index of the fragment defined",&
         usage="<INTEGER>",  type_of_var=integer_t, n_var=1)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
         description="Starting and ending atomic index defining one fragment must be provided",&
         usage="{Integer} {Integer}", type_of_var=integer_t, n_var=2, repeats=.FALSE.)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="MAP",&
         description="Provides the index of the fragment of the MIXED force_eval mapped on the"//&
         " locally defined fragment.",&
         usage="MAP <INTEGER>", type_of_var=integer_t, n_var=1, repeats=.FALSE.)
    CALL section_add_keyword(sub3section,keyword)
    CALL keyword_release(keyword)

    CALL section_add_subsection(sub2section,sub3section)
    CALL section_release(sub3section)
    CALL section_add_subsection(subsection,sub2section)
    CALL section_release(sub2section)

    CALL section_add_subsection(section,subsection)
    CALL section_release(subsection)

    CALL create_print_mix_section(subsection)
    CALL section_add_subsection(section,subsection)
    CALL section_release(subsection)
  END SUBROUTINE create_mix_section

! *****************************************************************************
!> \brief Create the print section for mixed
!> \param section the section to create
!> \author teo
! *****************************************************************************
  SUBROUTINE create_print_mix_section(section)
    TYPE(section_type), POINTER              :: section

    CHARACTER(len=*), PARAMETER :: routineN = 'create_print_mix_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.
    CPASSERT(.NOT.ASSOCIATED(section))
    CALL section_create(section,name="print",&
         description="Section of possible print options in MIXED env.",&
         n_keywords=0, n_subsections=1, repeats=.FALSE.)

    NULLIFY(print_key)

    CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
         description="Controls the printing of information during the evaluation of "//&
         "the mixed environment. ",&
         print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__")
    CALL section_add_subsection(section,print_key)
    CALL section_release(print_key)

    CALL cp_print_key_section_create(print_key,"DIPOLE",&
         description="Controls the printing of dipole information. "//&
         "Requires the DIPOLE calculation be active for all subforce_eval.", &
         print_level=medium_print_level,filename="__STD_OUT__")
    CALL section_add_subsection(section,print_key)
    CALL section_release(print_key)
  END SUBROUTINE create_print_mix_section

END MODULE input_cp2k_mixed
