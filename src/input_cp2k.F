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
  USE input_val_types, ONLY: val_type, no_t, logical_t, integer_t, real_t, char_t,&
       enum_t
  USE cp_output_handling, only: silent_print_level, low_print_level,&
       medium_print_level, high_print_level, debug_print_level
  USE xc_derivatives, only: create_xc_fun_section
  USE cp_control_types, only: xc_rho_no_smooth, xc_rho_nn50, xc_rho_nn10, &
         xc_rho_spline2_smooth,xc_rho_spline3_smooth,xc_deriv_pw,&
         xc_deriv_spline3, xc_deriv_spline2,&
         xc_deriv_nn50_smooth, xc_deriv_nn10_smooth, xc_deriv_spline2_smooth,&
         xc_deriv_spline3_smooth
  USE pw_spline_utils, only: no_precond,precond_spl3_aint, precond_spl3_1
  USE cp_output_handling, only: cp_print_key_section_create
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='input_cp2k'

  INTEGER, PARAMETER, public :: do_qs=1,do_fist=2,do_kg=3,do_qmmm=4,do_atom=5,&
         do_farming=6,do_test=8
  INTEGER, PARAMETER, PUBLIC :: nve_ensemble=101,nvt_ensemble=102,&
       npt_i_ensemble=103, npt_f_ensemble=104
  INTEGER, PARAMETER, PUBLIC :: md_init_start=101,md_pos_restart=102,&
       md_pov_vel_restart=103
  INTEGER, PARAMETER, public :: energy_run=1,energy_force_run=2,&
       md_run=3, geo_opt_run=4,mc_run=5
  INTEGER, PARAMETER, public :: atomic_guess=1,restart_guess=2,random_guess=3,&
       core_guess=4
  INTEGER, PARAMETER, public :: pw_interp=1,linear_interp=2,&
         spline3_nopbc_interp=3,spline3_pbc_interp=4 ! to move somewhere else
  public :: create_cp2k_root_section
  !!***
  !****************************************************************************
contains

  !!****f* input_cp2k/create_cp2k_root_section *
  !!
  !!   NAME
  !!     create_cp2k_root_section
  !!
  !!   FUNCTION
  !!     creates the input structure of the file used by cp2k
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - root_section: the input structure to be created
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE create_cp2k_root_section(root_section,error)
    TYPE(section_type), pointer :: root_section
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_cp2k_input_structure',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: section
    TYPE(keyword_type), pointer :: keyword

    failure=.false.

    CPPrecondition(.NOT.ASSOCIATED(root_section),cp_failure_level,routineP,error,failure)
    IF (.NOT.  failure) THEN
       CALL section_create(root_section,name="__ROOT__",&
            description="input file of cp2k",n_keywords=0, n_subsections=10,&
            repeats=.FALSE., required=.FALSE., error=error)
       nullify(section)

!       CALL create_pippo_section(section,error=error)
!       CALL section_add_subsection(root_section,section,error=error)
!       CALL section_release(section,error=error)

     CALL create_global_section(section,error=error)
     CALL section_add_subsection(root_section,section,error=error)
     CALL section_release(section,error=error)

     CALL create_md_section(section,error=error)
     CALL section_add_subsection(root_section,section,error=error)
     CALL section_release(section,error=error)

     CALL create_force_eval_section(section,error=error)
     CALL section_add_subsection(root_section,section,error=error)
     CALL section_release(section,error=error)

       nullify(keyword)
       CALL keyword_create(keyword,name="_DEFAULT_KEYWORD_",&
            description="ignores lines outside sections",&
            usage="", type_of_var=no_t,error=error)
       CALL section_add_keyword(root_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF
  END SUBROUTINE create_cp2k_root_section
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

  !!****f* input_cp2k/create_global_section [1.0] *
  !!
  !!   NAME
  !!     
  !!
  !!   FUNCTION
  !!     section to hold global settings for the whole program
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - section: the section to be created
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE create_global_section(section,error)
    TYPE(section_type), pointer :: section
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_global_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword

    failure=.false.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="global",&
         description="Section with general information regarding which kind "//&
         "of simulation to perform an parameters for the whole PROGRAM",&
         n_keywords=7, n_subsections=0, repeats=.FALSE., required=.true.,&
         error=error) ! to do: required=.false.

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
         usage="FFT_LIBRARY FFTSG|FFTW|FFTESSL...",&
         default_c_val="FFTSG",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="PRINT_LEVEL",&
         variants=(/"IOLEVEL"/),&
         description="How much output is written out.",&
         usage="PRINT_LEVEL HIGH",&
         default_i_val=medium_print_level,enum_c_vals=&
         (/"SILENT","LOW   ","MEDIUM","HIGH  ","DEBUG "/),&
         enum_i_vals=(/silent_print_level,low_print_level,medium_print_level,&
         high_print_level,debug_print_level/),error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="PROGRAM_NAME",&
         variants=(/"PROGRAM"/),&
         description="Which program should be run to compute forces,...",&
         usage="PROGRAM_NAME QUICKSTEP|FIST|KG|QMMM|...",&
         enum_c_vals=(/&
         "QUICKSTEP   ",&
         "QS          ",&
         "FIST        ",&
         "KIM_GORDON  ",&
         "KG          ",&
         "QMMM        ",&
         "ATOM        ",&
         "FARMING     ",&
         "TEST        " /),&
         enum_i_vals=(/ do_qs,do_qs,do_fist,do_kg,do_kg,do_qmmm,do_atom,&
         do_farming,do_test/),&
         default_i_val=do_qs, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="PROJECT_NAME",&
         variants=(/"PROJECT"/),&
         description="Name of the project (used to build the name of the "//&
         "trajectory, and other files generated by the program)",&
         usage="PROJECT_NAME myFancyProject",&
         default_c_val="out",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="RUN_TYPE",&
         description="Type of run that you want to perform Geometry "//&
         "optimization, md, montecarlo,...",&
         usage="RUN_TYPE MD|GEO_OPT|MC",&
         default_i_val=energy_force_run,&
         enum_c_vals=(/&
         "ENERGY               ",&
         "ENERGY_FORCE         ",&
         "MD                   ",&
         "MOLECULAR_DYNAMICS   ",&
         "GEO_OPT              ",&
         "GEOMETRY_OPTIMIZATION",&
         "MC                   ",&
         "MONTECARLO           "/),&
         enum_i_vals=(/energy_run,energy_force_run,md_run,md_run,&
         geo_opt_run,geo_opt_run,mc_run,mc_run/),&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="TRACE",&
         description="If a debug trace of the execution of the program should be written ",&
         usage="TRACE",&
         default_l_val=.FALSE.,lone_keyword_l_val=.true.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

  END SUBROUTINE create_global_section
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
    TYPE(section_type), pointer :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="dft",&
         description="parameter needed by dft programs",&
         n_keywords=3, n_subsections=4, repeats=.FALSE., required=.TRUE.,&
         error=error)

    nullify(keyword)
    CALL keyword_create(keyword, name="BASIS_SET_FILE_NAME",&
         description="The path of the basis set file",&
         usage="BASIS_SET_FILE_NAME ../BASISSET",&
         default_c_val="BASISSET",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="POTENTIAL_FILE_NAME",&
         description="The path to the pseudo potential file",&
         usage="POTENTIAL_FILE_NAME ../POTENTIAL",&
         default_c_val="POTENTIAL",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="RESTART_FILE_NAME",&
         description="The path to the pseudo potential file",&
         usage="RESTART_FILE_NAME ../RESTART",&
         default_c_val="RESTART",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    nullify(subsection)
    CALL create_scf_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

!FM    CALL create_qs_section(subsection,error)
!FM    CALL section_add_subsection(section, subsection, error=error)
!FM    CALL section_release(subsection,error=error)
    
!FM    CALL create_tddfpt_section(subsection,error)
!FM    CALL section_add_subsection(section, subsection, error=error)
!FM    CALL section_release(subsection,error=error)

    CALL create_mgrid_section(subsection,error=error)
    CALL section_add_subsection(section, subsection,error=error)
    CALL section_release(subsection,error=error)

    CALL create_xc_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

  END SUBROUTINE create_dft_section
  !***************************************************************************

!!****f* input_cp2k/create_xc_section [1.0] *
!!
!!   NAME
!!     create_xc_section
!!
!!   FUNCTION
!!     creates the input section for the xc part
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE create_xc_section(section,error)
  TYPE(section_type), pointer :: section
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='create_xc_section',&
        routineP=moduleN//':'//routineN
  TYPE(section_type), pointer :: subsection
  TYPE(keyword_type), pointer :: keyword

  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL section_create(section,"xc",&
         description="parameters needed calculate the xc potential",&
         n_keywords=5, n_subsections=2, repeats=.FALSE., required=.TRUE.,&
         error=error)

    NULLIFY(subsection,keyword)
    
    CALL keyword_create(keyword, name="density_cutoff",&
         description="The cutoff on the density used by the xc calculation",&
         usage="density_cutoff 1.e-11", default_r_val=1.0e-10_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="gradient_cutoff",&
         description="The cutoff on the gradient of the density used by the "//&
         "xc calculation",&
         usage="gradient_cutoff 1.e-11", default_r_val=1.0e-10_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="tau_cutoff",&
         description="The cutoff on tau used by the xc calculation",&
         usage="tau_cutoff 1.e-11", default_r_val=1.0e-10_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL section_create(subsection, name="grid",&
         description="The xc parameters used when calculating the xc on the grid",&
         n_keywords=5, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
         error=error)

    CALL keyword_create(keyword, name="xc_smooth_rho",&
         description="The density smoothing used for the xc calculation",&
         usage="xc_smooth_rho nn10", default_i_val=xc_rho_no_smooth,&
         enum_c_vals=(/ &
         "NONE   ",&
         "NN50   ",&
         "NN10   ",&
         "SPLINE2",&
         "NN6    ",&
         "SPLINE3",&
         "NN4    " /),&
         enum_i_vals=(/ xc_rho_no_smooth, xc_rho_nn50, xc_rho_nn10, &
         xc_rho_spline2_smooth,xc_rho_spline2_smooth,xc_rho_spline3_smooth,&
         xc_rho_spline3_smooth/),&
         error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="xc_deriv",&
         description="The method used to compute the derivatives",&
         usage="xc_deriv NN10_SMOOTH", default_i_val=xc_deriv_spline2,&
         enum_c_vals=(/&
         "PW            ",&
         "SPLINE3       ",&
         "SPLINE2       ",&
         "NN50_SMOOTH   ",&
         "NN10_SMOOTH   ",&
         "SPLINE2_SMOOTH",&
         "NN6_SMOOTH    ",&
         "SPLINE3_SMOOTH",&
         "NN4_SMOOTH    "/),&
         enum_i_vals=(/xc_deriv_pw, xc_deriv_spline3, xc_deriv_spline2,&
         xc_deriv_nn50_smooth, xc_deriv_nn10_smooth, xc_deriv_spline2_smooth,&
         xc_deriv_spline2_smooth, xc_deriv_spline3_smooth, xc_deriv_spline3_smooth/),&
         error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)
    
    CALL create_xc_fun_section(subsection,error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

 END IF
END SUBROUTINE create_xc_section
!***************************************************************************

!!****f* input_cp2k/create_scf_section *
!!
!!   NAME
!!     create_scf_section
!!
!!   FUNCTION
!!     creates the structure of the section with the dft scf parameters
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: will contain the scf section
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE create_scf_section(section,error)
  TYPE(section_type), pointer :: section
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='create_scf_section',&
        routineP=moduleN//':'//routineN
  TYPE(keyword_type), pointer :: keyword

  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL section_create(section,"scf",&
         description="parameters needed perform an scf run",&
         n_keywords=25, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
         error=error)

    nullify(keyword)
    CALL keyword_create(keyword, name="WRITE_RESTART_EACH",&
         description="How often the wavefunction restart file should be "//&
         "written during an scf (0=never).",usage="WRITE_RESTART_EACH 5",&
         default_i_val=10,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MAX_ITER_LUMO",&
         variants=(/"MAX_ITER_LUMOS"/),&
         description="The maximum number of iteration for the lumo computation",&
         usage="MAX_ITER_LUMO 100", default_i_val=2999,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="EPS_LUMO",&
         variants=(/"EPS_LUMOS"/),&
         description="target accuracy of the computation of the lumo energy",&
         usage="EPS_LUMO 1.e-6", default_r_val=1.0e-5_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="WORK_SYEVX",&
         description="limit the memory usage of SYEVX to a fraction of the maximum size",&
         usage="work_syevx 0.8", default_r_val=1.0_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="SMEAR",&
         description="Use a finite temperature like smearing of the occupation numbers",&
         usage="SMEAR 0.1", default_r_val=0._dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MAX_SCF",&
         description="Maximum number of SCF iteration to be performed for one optimization",&
         usage="MAX_SCF 200", default_i_val=50,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MAX_DIIS",&
         variants=(/"MAX_DIIS_BUFFER_SIZE"/),&
         description="Maximum number of DIIS vectors to be used",&
         usage="MAX_DIIS 3", default_i_val=4,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="LEVEL_SHIFT",&
         variants=(/"LSHIFT"/),&
         description="Use level shifting to improve convergence",&
         usage="LEVEL_SHIFT 0.1", default_r_val=0._dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="EPS_SCF",&
         description="target accuracy for the scf convergence",&
         usage="EPS_SCF 1.e-6", default_r_val=1.e-5_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="CHOLESKY",&
         description="If the cholesky method should be used for computing "//&
         "the inverse of S",&
         usage="CHOLESKY OFF", default_l_val=.true.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="JACOBI_THRESHOLD",&
         description="XXXXXXXX switch to jacobi if ...",&
         usage="JACOBI_THRESHOLD 1.0e-6", default_r_val=1.0E-7_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="EPS_JACOBI",&
         description="XXXXXXXX target jacobi precision",&
         usage="EPS_JACOBI 1.e-5", default_r_val=0._dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="EPS_EIGVAL",&
         description="Throw away linear combinations of basis functions with a small eigenvalue in S",&
         usage="EPS_EIGVAL 1.0", default_r_val=1.0e-5_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="EPS_DIIS",&
         description="Threshold on the convergence to start using DIAG/DIIS",&
         usage="EPS_DIIS 5.0e-2", default_r_val=0.1_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MIXING",&
         variants=(/"DENSITY_MIXING"/),&
         description="Fraction of new density to be mixed in",&
         usage="MIXING 0.1", default_r_val=0.4_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="SCF_GUESS",&
         description="Change the initial guess for the wavefunction.",&
         usage="SCF_GUESS RESTART", default_i_val=atomic_guess,&
         enum_c_vals=(/ "ATOMIC ", "RESTART", "RANDOM ", "CORE   "/),&
         enum_i_vals=(/atomic_guess,restart_guess,random_guess,core_guess/),&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="OT",&
         description="Use the orbital transformation method",&
         usage="OT", default_l_val=.FALSE.,lone_keyword_l_val=.true.,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="NROW_BLOCK",&
         description="sets the number of rows in a scalapack block",&
         usage="NROW_BLOCK 31", default_i_val=32,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="NCOL_BLOCK",&
         description="Sets the number of columns in a scalapack block",&
         usage="NCOL_BLOCK 31", default_i_val=32,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="ADDED_MOS",&
         description="Number of additional MOS added for each spin",&
         usage="ADDED_MOS", default_i_val=1,n_var=-1,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

 END IF
END SUBROUTINE create_scf_section
!***************************************************************************

  !!****f* input_cp2k/create_force_eval_section *
  !!
  !!   NAME
  !!     create_force_eval_section
  !!
  !!   FUNCTION
  !!     creates the force_eval section
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
  SUBROUTINE create_force_eval_section(section,error)
    TYPE(section_type), pointer :: section
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_force_eval_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword
    TYPE(section_type), pointer :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="force_eval",&
         description="parameters needed to calculate energy and forces and"//&
         " describe the system you want to analyze.",&
         n_keywords=0, n_subsections=5, repeats=.FALSE., required=.TRUE.,&
         error=error)


    nullify(subsection)
    CALL create_dft_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

!FM    CALL create_mm_section(subsection,error)
!FM    CALL section_add_subsection(section, subsection, error=error)
!FM    CALL section_release(subsection,error=error)
!FM    
!FM    CALL create_qmmm_section(subsection,error)
!FM    CALL section_add_subsection(section, subsection, error=error)
!FM    CALL section_release(subsection,error=error)
!FM
!FM    CALL create_metadyn_section(subsection,error)
!FM    CALL section_add_subsection(section, subsection, error=error)
!FM    CALL section_release(subsection,error=error)

    CALL create_subsys_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

  END SUBROUTINE create_force_eval_section
  !***************************************************************************

!!****f* input_cp2k/create_md_section *
!!
!!   NAME
!!     create_md_section
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: will contain the md section 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE create_md_section(section,error)
  TYPE(section_type), pointer :: section
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='create_md_section',&
        routineP=moduleN//':'//routineN
  TYPE(keyword_type), pointer :: keyword

  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  if (.not. failure) then
    CALL section_create(section,name="md_new",&
         description="parameters needed perform an md run",&
         n_keywords=5, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
         error=error)

    nullify(keyword)
    CALL keyword_create(keyword, name="ensemble",&
         description="The ensemble that you want to sample with the md."//&
         " npt_i uses an isotropic resize of the cell, npt_f a flexible cell.",&
         usage="ensemble nve",&
         default_i_val=nve_ensemble,&
         enum_c_vals=(/"nve  ","nvt  ","npt_i","npt_f"/),&
         enum_i_vals=(/nve_ensemble,nvt_ensemble,npt_i_ensemble,&
         npt_f_ensemble/),error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="steps",&
         description="The number of steps to be performed",&
         usage="steps 100", default_i_val=3,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="timestep",&
         description="The length of an integration step (in fs).",&
         usage="timestep 1.0", default_r_val=0.5_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="restart",&
         description="The way the trajectory should be done (init: from the "//&
         "initial (input file) positions, with velocities from a boltzmann "//&
         "distribution; pos: using the latest positions stored in the restart"//&
         " file, with with velocities from a boltzmann distribution"//&
         "; posvel: using both the latest positions and velocites stored in "//&
         "the restart file", usage="RESTART posvel", &
         default_i_val=md_init_start,&
         enum_c_vals=(/&
         "init  ",&
         "pos   ",&
         "posvel"/),&
         enum_i_vals=(/md_init_start,md_pos_restart, md_pov_vel_restart/),&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    
    CALL keyword_create(keyword, name="temperature",&
         description="The temperature in K used to initialize the velocities"//&
         "with init and pos restart, and in the npt/nvt simulations", &
         usage="temperature 325.0",&
         default_r_val=300.0_dp, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    
    CALL keyword_create(keyword, name="temp_tol",&
         variants=(/"temp_to              ", "temperature_tolerance"/),&
         description="The maximum accepted oscillations of the temperature"//&
         "from the wanted temperature before rescaling the velocites."//&
         "If it is 0 no rescaling is performed", &
         usage="temp_tol 0.0", default_r_val=0.0_dp, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    
  end if
END SUBROUTINE create_md_section
!***************************************************************************

!!****f* input_cp2k/create_subsys_section [1.0] *
!!
!!   NAME
!!     create_subsys_section
!!
!!   FUNCTION
!!     creates the structure of a subsys, i.e. a full set of 
!!     atoms+mol+bounds+cell
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE create_subsys_section(section,error)
  TYPE(section_type), pointer :: section
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='create_subsys_section',&
        routineP=moduleN//':'//routineN
  TYPE(section_type), pointer :: subsection

  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  if (.not. failure) then
    CALL section_create(section,name="subsys",&
         description="a subsystem: coordinates, topology, molecules and cell",&
         n_keywords=5, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
         error=error)

    nullify(subsection)
!FM    CALL create_cell_section(subsection,error=error)
!FM    CALL section_add_subsection(section,subsection,error=error)
!FM    CALL section_release(subsection,error=error)
!FM
!FM    CALL create_coord_section(subsection,error=error)
!FM    CALL section_add_subsection(section,subsection,error=error)
!FM    CALL section_release(subsection,error=error)
!FM
!FM    CALL create_topology_section(subsection,error=error)
!FM    CALL section_add_subsection(section,subsection,error=error)
!FM    CALL section_release(subsection,error=error)
!FM
!FM    CALL create_kind_section(subsection,error=error)
!FM    CALL section_add_subsection(section,subsection,error=error)
!FM    CALL section_release(subsection,error=error)
     
  end if
END SUBROUTINE create_subsys_section
!***************************************************************************


!!****f* input_cp2k/create_mgrid_section *
!!
!!   NAME
!!     create_mgrid_section
!!
!!   FUNCTION
!!     creates the multigrid
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE create_mgrid_section(section,error)
  TYPE(section_type), pointer :: section
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='create_mgrid_section',&
        routineP=moduleN//':'//routineN
  TYPE(section_type), POINTER :: subsection,print_key
  TYPE(keyword_type), pointer :: keyword
  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
    CALL section_create(section,name="mgrid",&
         description="multigrid information",&
         n_keywords=5, n_subsections=1, repeats=.FALSE., required=.false.,&
         error=error)
    nullify(keyword)
    CALL keyword_create(keyword, name="ngrids",&
         description="The number of multigrids to use",&
         usage="ngrids 1", default_i_val=100, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="cutoff",&
         description="The cutoff (from the highest grid to the lowest)",&
         usage="cutoff 300", default_r_val=280._dp, n_var=-1, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="progression_factor",&
         description="Factor used to find the cutoff of the multigrids that"//&
         " where not given explicetly",&
         usage="progression_factor 4", default_i_val=3, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="commensurate",&
         description="If the grids should be commensurate. If true overrides "//&
         "the progression factor and the cutoffs of the sub grids",&
         usage="commensurate", default_l_val=.FALSE.,&
              lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    nullify(subsection)
    CALL section_create(subsection,name="interpolator",&
         description="kind of interpolation used between the multigrids",&
         n_keywords=5, n_subsections=0, repeats=.FALSE., required=.false.,&
         error=error)
    CALL keyword_create(keyword, name="kind",&
         description="the interpolator to use",&
         usage="kind spline3",&
         default_i_val=pw_interp,&
         enum_c_vals=(/&
         "pw           ",&
         "linear       ",&
         "spline3_nopbc",&
         "spline3      "/),&
         enum_i_vals=(/pw_interp,linear_interp,&
         spline3_nopbc_interp,spline3_pbc_interp/),error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="aint_precond",&
         description="the approximate inverse to use to get the starting point"//&
         " for the linear solver of the spline3 methods",&
         usage="kind spline3",&
         default_i_val=pw_interp,&
         enum_c_vals=(/&
         "copy               ",&
         "spl3_nopbc_aint    ",&
         "spl3_nopbc_precond1"/),&
         enum_i_vals=(/no_precond,precond_spl3_aint, precond_spl3_1/),&
         error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="precond",&
         description="The preconditioner used"//&
         " for the linear solver of the spline3 methods",&
         usage="kind spline3",&
         default_i_val=pw_interp,&
         enum_c_vals=(/&
         "copy               ",&
         "spl3_nopbc_aint    ",&
         "spl3_nopbc_precond1"/),&
         enum_i_vals=(/no_precond,precond_spl3_aint, precond_spl3_1/),&
         error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="eps_x",&
         description="accuracy on the solution for spline3 the interpolators",&
         usage="eps_x 1.e-15", default_r_val=1.e-10_dp, error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="eps_r",&
         description="accuracy on the residual for spline3 the interpolators",&
         usage="eps_r 1.e-15", default_r_val=1.e-10_dp, error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="max_iter",&
         description="the maximum number of iterations",&
         usage="max_iter 200", default_i_val=100, error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    nullify(print_key)
    CALL cp_print_key_section_create(print_key,"conv_info",&
         description="if convergence information about the linear solver"//&
         " of the spline methods should be printed", &
         print_level=medium_print_level,each=(/10/),filename="__STD_OUT__",&
         error=error)
    CALL section_add_subsection(subsection,print_key,error=error)
    CALL section_release(print_key,error=error)
    CALL cp_print_key_section_create(print_key,"spl_coeffs",&
         description="outputs a cube with the coefficents calculated for "//&
         "the spline interpolation", &
         print_level=debug_print_level,each=(/1/),&
         error=error)
    CALL section_add_subsection(subsection,print_key,error=error)
    CALL section_release(print_key,error=error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)
 END IF
END SUBROUTINE create_mgrid_section
!***************************************************************************

END MODULE input_cp2k
