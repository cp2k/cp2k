!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/input_cp2k [1.0] *
!!
!!   NAME
!!     input_cp2k
!!
!!   FUNCTION
!!     builds the input structure for cp2k
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE input_cp2k
  use cp_log_handling, only: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log
  use cp_error_handling, only: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger, cp_a_l, cp_error_check
  use kinds, only: dp
  use cp_para_types, only: cp_para_env_type
  use timings, only: timeset, timestop
  USE input_section_types, ONLY: section_type, section_create, &
       section_add_keyword, section_add_subsection, section_release
  USE input_keyword_types, ONLY: keyword_type, keyword_create,&
       keyword_release
  USE input_file_types, ONLY: input_structure_type, input_structure_create
  USE input_val_types, ONLY: val_type, no_t, logical_t, integer_t, real_t, char_t,&
       enum_t
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='input_cp2k'

  public :: create_cp2k_input_structure
  !!***
  !****************************************************************************
contains

  !!****f* input_cp2k/create_cp2k_input_structure *
  !!
  !!   NAME
  !!     create_cp2k_input_structure
  !!
  !!   FUNCTION
  !!     creates the input structure of the file used by cp2k
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - input_structure: the input structure to be created
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE create_cp2k_input_structure(input_structure,error)
    TYPE(input_structure_type), pointer :: input_structure
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_cp2k_input_structure',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: root_section,section
    TYPE(keyword_type), pointer :: keyword

    failure=.false.

    CPPrecondition(.NOT.ASSOCIATED(input_structure),cp_failure_level,routineP,error,failure)
    CALL input_structure_create(input_structure,&
         description="input file of cp2k",&
         n_sections=10,error=error)
    root_section => input_structure%root_section
    IF (.NOT.  failure) THEN
       nullify(section)

       CALL create_pippo_section(section,error=error)
       CALL section_add_subsection(root_section,section,error=error)
       CALL section_release(section,error=error)

!FM     CALL create_cp2k_section(section,error=error)
!FM     CALL section_add_subsection(root_section,section,error=error)
!FM     CALL section_release(section,error=error)
!FM
!FM     CALL create_dft_section(section,error=error)
!FM     CALL section_add_subsection(root_section,section,error=error)
!FM     CALL section_release(section,error=error)

       nullify(keyword)
       CALL keyword_create(keyword,name="_DEFAULT_KEYWORD_",&
            description="ignores lines outside sections",&
            usage="", type_of_var=no_t,error=error)
       CALL section_add_keyword(root_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF
  END SUBROUTINE create_cp2k_input_structure
  !***************************************************************************

  !!****f* input_cp2k/create_pippo_section *
  !!
  !!   NAME
  !!     create_pippo_section
  !!
  !!   FUNCTION
  !!     creates a test input section with many combinations of inputs
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - section: the section that is created
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE create_pippo_section(section,error)
    TYPE(section_type), POINTER :: section
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_pippo_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: subsection
    TYPE(keyword_type), pointer :: keyword

    failure=.false.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT.  failure) THEN
       nullify(section)
       CALL section_create(section,name="pippo",&
            description="Test section"//&
            " bla bla bla bla bla bla blah bla blah bla bla bla bla bla bla blah bla blah",&
            n_keywords=7, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       nullify(keyword)
       CALL keyword_create(keyword, name="SINGLE_logical",&
            description="a logical value ", usage="SINGLE_LOGICAL yes",&
            default_l_val=.FALSE.,lone_keyword_l_val=.true.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="2_logicals",&
            description="2 logical values ", usage="2_LOGICALS yes F",&
            n_var=2,type_of_var=logical_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="n_logicals",&
            description="many logical values ", usage="n_LOGICALS yes F T",&
            n_var=-1,type_of_var=logical_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="r_logical",&
            description="repeated logical value ", usage="r_LOGICAL yes",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,repeats=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       nullify(subsection)
       CALL section_create(subsection,name="int",&
            description="Integer test subsection"//&
            " bla bla bla bla bla bla blah bla blah bla bla bla bla bla bla blah bla blah",&
            n_keywords=7, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="SINGLE_integer",&
            description="a integer value ", usage="SINGLE_INTEGER yes",&
            default_i_val=0,lone_keyword_i_val=1,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="2_integers",&
            description="2 integer values ", usage="2_INTEGERS 3 4",&
            n_var=2,type_of_var=integer_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="n_integers",&
            description="many integer values ", usage="n_INTEGERS 1 3 6",&
            n_var=-1,type_of_var=integer_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="r_integer",&
            description="repeated integer value ", usage="r_INTEGER 4",&
            default_i_val=1,lone_keyword_i_val=3,repeats=.TRUE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       nullify(subsection)
       CALL section_create(subsection,name="real",&
            description="Reals test subsection"//&
            " bla bla bla bla bla bla blah bla blah bla bla bla bla bla bla blah bla blah",&
            n_keywords=7, n_subsections=0, repeats=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="SINGLE_real",&
            description="a real value ", usage="SINGLE_REAL yes",&
            default_r_val=0._dp,lone_keyword_r_val=1._dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="2_reals",&
            description="2 real values ", usage="2_REALS 3 4",&
            n_var=2,type_of_var=real_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="n_reals",&
            description="many real values ", usage="n_REALS 1 3 6",&
            n_var=-1,type_of_var=real_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="r_real",&
            description="repeated real value ", usage="r_REAL 4",&
            default_r_val=1._dp,lone_keyword_r_val=3._dp,repeats=.TRUE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       nullify(subsection)
       CALL section_create(subsection,name="char",&
            description="Chars test subsection"//&
            " bla bla bla bla bla bla blah bla blah bla bla bla bla bla bla blah bla blah",&
            n_keywords=7, n_subsections=0, repeats=.FALSE., required=.false.,&
            error=error)

       CALL keyword_create(keyword, name="SINGLE_char",&
            description="a char value ", usage="SINGLE_CHAR a_word",&
            default_c_val="dv",lone_keyword_c_val="lv",error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="2_chars",&
            description="2 char values ", usage='2_CHARS p "bla bla"',&
            n_var=2,type_of_var=char_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="n_chars",&
            description="many char values ", usage='n_CHARS t 1 tanti "t r"',&
            n_var=-1,type_of_var=char_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="r_char",&
            description="repeated char value ", usage="r_CHAR 4",&
            default_c_val="default",lone_keyword_c_val="value0",repeats=.TRUE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       nullify(subsection)
       CALL section_create(subsection,name="enum",&
            description="Enum test subsection"//&
            " bla bla bla bla bla bla blah bla blah bla bla bla bla bla bla blah bla blah",&
            n_keywords=7, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="SINGLE_enum",&
            description="a enum value ", usage="SINGLE_ENUM pippo0",&
            default_i_val=0,lone_keyword_i_val=1,&
            enum_c_vals=(/"pippo0","pippo1","p1    ","p0    "/),&
            enum_i_vals=(/ 0,1,1,0 /),error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="2_enums",&
            description="2 enum values ", usage="2_ENUMS p1 p0",&
            n_var=2,type_of_var=enum_t,&
            enum_c_vals=(/"pippo0","pippo1","p1    ","p0    "/),&
            enum_i_vals=(/ 0,1,1,0 /),error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="n_enums",&
            description="many enum values ", usage="n_ENUMS pippo0 pippo1",&
            enum_c_vals=(/"pippo0","pippo1","p1    ","p0    "/),&
            enum_i_vals=(/ 0,1,1,0 /),n_var=-1,type_of_var=enum_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="r_enum",&
            description="repeated enum value ", usage="r_ENUM p1",&
            enum_c_vals=(/"pippo0","pippo1","p3    ","p0    "/),&
            enum_i_vals=(/ 0,1,3,0 /),enum_strict=.FALSE.,&
            default_i_val=1,lone_keyword_i_val=3,&
            repeats=.TRUE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_pippo_section
  !***************************************************************************

  !!****f* input_cp2k/create_cp2k_section [1.0] *
  !!
  !!   NAME
  !!     
  !!
  !!   FUNCTION
  !!     adds the generic sections: CP2K, cell
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE create_cp2k_section(section,error)
    TYPE(section_type), pointer :: section
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_cp2k_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword

    failure=.false.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="cp2k",&
         description="Section with general information regarding which kind "//&
         "of simulation to perform an parameters for the whole PROGRAM",&
         n_keywords=7, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
         error=error)

    nullify(keyword)
    CALL keyword_create(keyword, name="BLACS_GRID",&
         description="how to distribute the processors on the 2d grid needed "//&
         "by BLACS (and thus SCALAPACK)", usage="BLACS_GRID SQUARE|COL|COLUMN|ROW",&
         default_c_val="BLACS_GRID_SQUARE",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="FFT_LIBRARY",&
         variants=(/"FFT_LIB","FFTLIB "/),&
         description="which fft library to use to perform the ffts",&
         usage="FFT_LIBRARY FFTSG|FFTW|...",&
         default_c_val="FFTSG",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="IOLEVEL",&
         variants=(/"PRINT_LEVEL"/),&
         description="How much output is written out: 0: few, 4: average, 7: more, 10: lots",&
         usage="IOLEVEL 3",&
         default_i_val=0,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="PROGRAM_NAME",&
         variants=(/"PROGRAM"/),&
         description="Which program should be run to compute forces,...",&
         usage="PROGRAM_NAME QUICKSTEP|FIST|KG|QMMM|...",&
         default_c_val="QUICKSTEP",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="PROJECT_NAME",&
         variants=(/"PROJECT"/),&
         description="Name of the project (used to build the name of the "//&
         "trajectory, and other files generated by the program)",&
         usage="PROJECT_NAME myFancyProject",&
         default_c_val="",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="RUN_TYPE",&
         description="Type of run that you want to perform Geometry "//&
         "optimization, md, montecarlo,...",&
         usage="RUN_TYPE MD|GEO_OPT|MC",&
         default_c_val="",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="TRACE",&
         description="If a debug trace of the execution of the program should be written ",&
         usage="TRACE",&
         default_l_val=.FALSE.,lone_keyword_l_val=.true.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

  END SUBROUTINE create_cp2k_section
  !***************************************************************************

  !!****f* input_cp2k/create_dft_section *
  !!
  !!   NAME
  !!     create_dft_section
  !!
  !!   FUNCTION
  !!     creates the dft section
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - section: the section to be created
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE create_dft_section(section,error)
    TYPE(section_type), pointer :: section
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_dft_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="dft",&
         description="parameter needed by any dft program for the xc functional",&
         n_keywords=16, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
         error=error)

    nullify(keyword)
    CALL keyword_create(keyword, name="functional",&
         description="the functional to be used",&
         usage="FUNCTIONAL BLYP|PADE|...",&
         default_c_val="PADE",repeats=.true.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="density_cutoff",&
         description="The valoue below witch the density is ignored for xc evaluation",&
         usage="DENSITY_CUTOFF 1.e-10",&
         default_r_val=1.e-10_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

  END SUBROUTINE create_dft_section
  !***************************************************************************

END MODULE input_cp2k
