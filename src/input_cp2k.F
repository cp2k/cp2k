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
  USE cp_log_handling, ONLY: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log
  USE cp_error_handling, ONLY: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger, cp_a_l, cp_error_check
  USE kinds, ONLY: dp, dbl
  USE cp_para_types, ONLY: cp_para_env_type
  USE timings, ONLY: timeset, timestop
  USE input_section_types, ONLY: section_type, section_create, &
       section_add_keyword, section_add_subsection, section_release,&
       section_describe_html, section_describe_index_html, section_vals_type,&
       section_vals_get_subs_vals, section_vals_get
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
       xc_deriv_spline3_smooth,sic_none, sic_mauri_us, sic_ad, sic_mauri_spz,&
       no_excitations,tddfpt_excitations
  USE pw_spline_utils, only: no_precond,precond_spl3_aint, precond_spl3_1
  USE cp_output_handling, only: cp_print_key_section_create
  USE greens_fn, ONLY: periodic3d,& 
       analytic2d,& 
       analytic1d,& 
       analytic0d,& 
       mt2d,&       
       mt1d,&       
       mt0d,&       
       hockney2d,&
       hockney1d,&
       hockney0d 
  USE string_utilities, ONLY: s2a
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
       md_run=3, geo_opt_run=4,mc_run=5,electronic_spectra_run=6
  INTEGER, PARAMETER, public :: atomic_guess=1,restart_guess=2,random_guess=3,&
       core_guess=4
  INTEGER, PARAMETER, public :: pw_interp=1,linear_interp=2,&
       spline3_nopbc_interp=3,spline3_pbc_interp=4 ! to move somewhere else


  integer, parameter, public :: do_gapw_gcs=1,&
       do_gapw_gct=2,&
       do_gapw_log=3
  integer, parameter, public :: do_pwgrid_spherical=1,&
       do_pwgrid_ns_fullspace=2,&
       do_pwgrid_ns_halfspace=3
  INTEGER, PARAMETER, PUBLIC :: do_etrxpl_use_guess=0,&
       do_etrxpl_use_prev_p=1,&
       do_etrxpl_use_prev_rho_r=2,&
       do_etrxpl_linear_wf=3,&
       do_etrxpl_linear_p=4,&
       do_etrxpl_linear_ps=5,&
       do_etrxpl_use_prev_wf=6,&
       do_etrxpl_ps=7
  INTEGER, PARAMETER, PUBLIC :: do_method_gapw=1,&
       do_method_gpw=2,&
       do_method_kg_pol=3,&
       do_method_kg_nopol=4,&
       do_method_kg_gpw=5
  integer, parameter, public :: do_tddfpt_singlet=0,&
                                do_tddfpt_triplet=1,&
                                do_tddfpt_davidson=0,&
                                do_tddfpt_lanczos=1,&
                                oe_none=0, oe_saop=1, oe_sic=2
  integer, parameter, public :: use_unit_angstrom=0,&
       use_unit_bohr=1,&
       use_unit_scaled_angstrom=2,&
       use_unit_scaled_bohr=3
  integer, parameter, public :: do_rsgrid_distributed=0,&
       do_rsgrid_replicated=1

  integer, parameter, public :: do_extrpl_ps=0,&
       do_extrpl_use_guess=1,&
       do_extrpl_use_prev_p=2,&
       do_extrpl_use_prev_wf=3,&
       do_extrpl_use_prev_rho_r=4,&
       do_extrpl_linear_wf=5,&
       do_extrpl_linear_p=6,&
       do_extrpl_linear_ps=7
  integer, parameter, public :: do_qmmm_none=0,&
                                do_qmmm_grid=1,&
                                do_qmmm_gauss=2

  integer, parameter, public :: do_ewald=0,&
                                do_pme=1,&
                                do_spme=2

  INTEGER, PARAMETER :: ot_mini_sd=1,ot_mini_cg=2,ot_mini_diis=3
  INTEGER, PARAMETER :: ot_ls_none=1,ot_ls_2pnt=2,ot_ls_3pnt=3,ot_ls_gold=4
  INTEGER, PARAMETER :: ot_precond_none=0, &
                        ot_precond_full_single=1, &
                        ot_precond_full_kinetic=2, &
                        ot_precond_s_inverse=3, &
                        ot_precond_sparse_diag=4


  PUBLIC :: create_cp2k_root_section, check_cp2k_input, parsed_cp2k_input
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

!!****f* input_cp2k/check_cp2k_input *
!!
!!   NAME
!!     check_cp2k_input
!!
!!   FUNCTION
!!     performs further checks on an input that parsed successfully
!!
!!   NOTES
!!     at the moment does nothing
!!
!!   ARGUMENTS
!!     - input_file: the parsed input
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE check_cp2k_input(input_file,error)
    TYPE(section_vals_type), pointer :: input_file
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='check_cp2k_input',&
         routineP=moduleN//':'//routineN

    failure=.false.

    CPPrecondition(ASSOCIATED(input_file),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPrecondition(input_file%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT. failure) THEN
! do the checks!!! (Yes, I agree, but, aren't they done elsewhere already ?)
    END IF
  END SUBROUTINE check_cp2k_input
!***************************************************************************

!!****f* input_cp2k/parsed_cp2k_input *
!!
!!   NAME
!!     parsed_cp2k_input
!!
!!   FUNCTION
!!     utility function to ease the transition to the new input.
!!     returns true if the new input was parsed
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - input_file: the parsed input file
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  FUNCTION parsed_cp2k_input(input_file,error) RESULT(res)
    TYPE(section_vals_type), pointer :: input_file
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    logical :: res

    logical :: failure
    character(len=*), parameter :: routineN='parsed_cp2k_input',&
         routineP=moduleN//':'//routineN
    TYPE(section_vals_type), POINTER ::glob_section

    failure=.false.

    res=ASSOCIATED(input_file)
    IF (res) THEN
       CPPrecondition(input_file%ref_count>0,cp_failure_level,routineP,error,failure)
       IF (.NOT.failure) THEN
          glob_section => section_vals_get_subs_vals(input_file,"GLOBAL",&
               error=error)
          CALL section_vals_get(glob_section,explicit=res,error=error)
       END IF
    END IF
  END FUNCTION parsed_cp2k_input
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
         "by BLACS (and thus SCALAPACK)", usage="BLACS_GRID SQUARE",&
         default_c_val="BLACS_GRID_SQUARE",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="FFT_LIBRARY",&
         variants=s2a("FFT_LIB","FFTLIB"),&
         description="which fft library to use to perform the ffts",&
         usage="FFT_LIBRARY FFTW",&
         default_c_val="FFTSG",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="PRINT_LEVEL",&
         variants=(/"IOLEVEL"/),&
         description="How much output is written out.",&
         usage="PRINT_LEVEL HIGH",&
         default_i_val=medium_print_level,enum_c_vals=&
          s2a("SILENT","LOW","MEDIUM","HIGH","DEBUG"),&
         enum_desc=s2a( "Almost no output",&
                        "Little output", "Quite some output", "Lots of output",&
                        "Everything is written out, useful for debugging purposes only"),&
         enum_i_vals=(/silent_print_level,low_print_level,medium_print_level,&
         high_print_level,debug_print_level/),error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="PROGRAM_NAME",&
         variants=(/"PROGRAM"/),&
         description="Which program should be run to compute forces,...",&
         usage="PROGRAM_NAME QUICKSTEP",&
         enum_c_vals=s2a("QUICKSTEP","QS","FIST","KIM_GORDON","KG",&
         "QMMM", "ATOM", "FARMING", "TEST"),&
         enum_desc=s2a("Electronic structure methods (DFT, ...)","Idem",&
                       "Molecular Mechanics","Kim Gordon models (DFT)",&
                       "Idem","Hybrid quantum classical","Atomic calculations",&
                       "Runs N independent jobs in a single run",&
                       "Do some benchmarking and testing"),&
         enum_i_vals=(/do_qs,do_qs,do_fist,do_kg,do_kg,do_qmmm,do_atom,&
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
         usage="RUN_TYPE MD",&
         default_i_val=energy_force_run,&
         enum_c_vals=s2a( "ENERGY", "ENERGY_FORCE", "MD", "MOLECULAR_DYNAMICS",&
         "GEO_OPT", "GEOMETRY_OPTIMIZATION", "MC", "MONTECARLO",&
         "SPECTRA", "ELECTRONIC_SPECTRA"),&
         enum_i_vals=(/energy_run,energy_force_run,md_run,md_run,&
         geo_opt_run,geo_opt_run,mc_run,mc_run,&
         electronic_spectra_run,electronic_spectra_run/),&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="ECHO_INPUT",&
         description="If the input should be echoed to the output with all the "//&
         "defaults made explicit",&
         usage="ECHO_INPUT NO",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
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
         description="Name of the basis set file, may include a path",&
         usage="BASIS_SET_FILE_NAME ../BASISSET",&
         default_c_val="BASISSET",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="POTENTIAL_FILE_NAME",&
         description="Name of the pseudo potential files, may include a path",&
         usage="POTENTIAL_FILE_NAME ../POTENTIAL",&
         default_c_val="POTENTIAL",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="LSD",&
         description="If the calculation should open shell",&
         usage="lsd",&
         default_l_val=.FALSE.,lone_keyword_l_val=.true.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="RESTRICTED",&
         description="If the calculation should be restricted",&
         usage="RESTRICTED",&
         default_l_val=.FALSE.,lone_keyword_l_val=.true.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="MULTIPLICITY",&
         description="Two times the total spin plus one.",&
         usage="MULTIPLICITY 1",&
         default_i_val=0,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="CHARGE",&
         description="The total charge of the system",&
         usage="CHARGE -1",&
         default_i_val=0,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="EXCITATIONS",&
         description="If excitations should be calculated",&
         usage="EXCITATIONS",&
         enum_c_vals=s2a("NONE","TDLR","TDDFPT"),&
         enum_i_vals=(/ no_excitations, tddfpt_excitations, &
         tddfpt_excitations/),&
         default_i_val=no_excitations, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    nullify(subsection)
    CALL create_scf_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_qs_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_tddfpt_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_mgrid_section(subsection,error=error)
    CALL section_add_subsection(section, subsection,error=error)
    CALL section_release(subsection,error=error)

    CALL create_xc_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_sic_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

  END SUBROUTINE create_dft_section
!***************************************************************************

!!****f* input_cp2k/create_qs_section [1.0] *
!!
!!   NAME
!!     create_qs_section
!!
!!   FUNCTION
!!     creates the input section for the qs part
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
!!     teo
!!
!!*** **********************************************************************
SUBROUTINE create_qs_section(section,error)
  ! Arguments
  TYPE(section_type), pointer :: section
  type(cp_error_type), optional, intent(inout) :: error
  ! Local Variables
  logical :: failure
  character(len=*), parameter :: routineN='create_qs_section',&
        routineP=moduleN//':'//routineN
  TYPE(keyword_type), pointer :: keyword

  failure=.false.

  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL section_create(section,"qs",&
          description="parameters needed to set up the Quickstep framework",&
          n_keywords=34, n_subsections=2, repeats=.FALSE., required=.FALSE.,&
          error=error)

     NULLIFY(keyword)
    
    CALL keyword_create(keyword, name="RESTART_FILE_NAME",&
         description="The path to the pseudo potential file",&
         usage="RESTART_FILE_NAME ../RESTART",&
         default_c_val="RESTART",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
     ! Reals
     CALL keyword_create(keyword, name="EPS_DEFAULT",&
          description="Try setting all EPS_xxx to values leading to an energy correct up to EPS_DEFAULT",&
          usage="EPS_DEFAULT real", default_r_val=1.0E-10_dp,&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)     

     CALL keyword_create(keyword, name="EPS_CORE_CHARGE",&
          description="Precision for mapping the core charges",&
          usage="EPS_CORE_CHARGE real", type_of_var=real_t, n_var=1,&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)     

     CALL keyword_create(keyword, name="EPS_GVG_RSPACE",&
          variants=(/"EPS_GVG"/),&
          description="Sets precision of the realspace KS matrix element integration",&
          usage="EPS_GVG_RSPACE real", type_of_var=real_t, n_var=1,&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="EPS_PGF_ORB",&
          description="Sets precision of the overlap matrix elements",&
          usage="EPS_PGF_ORB real", type_of_var=real_t, n_var=1,&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)


     CALL keyword_create(keyword, name="EPS_PPL",&
          description="Sets precision of the local part of the pseudo potential",&
          usage="EPS_PPL real", type_of_var=real_t, n_var=1,&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="EPS_PPNL",&
          description="Sets precision of the non-local part of the pseudo potential",&
          usage="EPS_PPNL real", type_of_var=real_t, n_var=1,&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="EPS_RHO",&
          description="Sets precision of the density mapping on the grids",&
          usage="EPS_RHO real", type_of_var=real_t, n_var=1,&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="EPS_RHO_RSPACE",&
          description="Sets precision of the density mapping in rspace",&
          usage="EPS_RHO_RSPACE real", type_of_var=real_t, n_var=1,&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="EPSFIT",&
          variants=(/"EPS_FIT"/),&
          description="GAPW : precision to give the extention of a hard gaussian ",&
          usage="EPSFIT real", default_r_val=1.0E-4_dp,&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="EPSISO",&
          variants=(/"EPS_ISO"/),&
          description="GAPW : precision to determine an isolated projector",&
          usage="EPSISO real", default_r_val=1.0E-12_dp,&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="EPSRHO0",&
          variants=s2a("EPSVRHO0","EPS_VRHO0"),&
          description="GAPW : precision to determine the range of V(rho0-rho0soft)",&
          usage="EPSRHO0 real", default_r_val=1.0E-6_dp,&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="ALPHA0_HARD",&
          variants=(/"ALPHA0_H"/),&
          description="GAPW : Exponent for hard compensation charge",&
          usage="ALPHA0_HARD real", default_r_val=0.0_dp,&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="ALPHA0_SOFT",&
          variants=(/"ALPHA0_S"/),&
          description="GAPW : Exponent for soft compensation charge",&
          usage="ALPHA0_SOFT real", default_r_val=0.0_dp,&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     ! Logicals
     CALL keyword_create(keyword, name="MAP_PAA",&
          description="Put diagonal density matrix elements on the highest cutoff multigrid",&
          usage="MAP_PAA",&
          default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="MAP_CONSISTENT",&
          description="Guarantee an exact derivative (Hks) of the energy",&
          usage="MAP_CONSISTENT",&
          default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="MULLIKEN_RESTRAINT",&
          description="Employ a restraint on the mulliken charges",&
          usage="MULLIKEN_RESTRAINT",&
          default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)     

     ! Integers
     CALL keyword_create(keyword, name="LMAXN1",&
          variants=(/"LMAXRHO1"/),&
          description="GAPW : max L number for espansion of the atomic densities in spherical gaussians",&
          usage="LMAXN1 integer",&
          default_i_val=6,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="LMAXN0",&
          variants=(/"LMAXRHO0"/),&
          description="GAPW : max L number for the expansion compensation densities in spherical gaussians",&
          usage="LMAXN0 integer",&
          default_i_val=2,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)    
     
     CALL keyword_create(keyword, name="LADDN0",&
          description="GAPW : XXXXXXXXXX",&
          usage="LADDN0 integer",&
          default_i_val=99,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     ! Characters
     CALL keyword_create(keyword, name="QUADRATURE",&
          description="GAPW: algorithm to construct the atomic radial grids",&
          usage="QUADRATURE GC_SIMPLE",&
          enum_c_vals=s2a("GC_SIMPLE","GC_TRANSFORMED","GC_LOG"),&
          enum_i_vals=(/ do_gapw_gcs,do_gapw_gct,do_gapw_log/),&
          default_i_val=do_gapw_log, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="PW_GRID",&
          description="What kind of PW_GRID should be employed",&
          usage="PW_GRID NS-FULLSPACE",&
          enum_c_vals=s2a("SPHERICAL","NS-FULLSPACE","NS-HALFSPACE"),&
          enum_desc=s2a("- not tested"," tested"," - not tested"),&
          enum_i_vals=(/ do_pwgrid_spherical, do_pwgrid_ns_fullspace,do_pwgrid_ns_halfspace/),&
          default_i_val=do_pwgrid_ns_fullspace, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="RS_GRID",&
          description="Parallellisation strategy of the realspace grids",&
          usage="RS_GRID DISTRIBUTED",&
          enum_c_vals=s2a("DISTRIBUTED","REPLICATED"),&
          enum_desc=s2a("should be faster in principle, is often not, rarely tested", "default"),&
          enum_i_vals=(/ do_rsgrid_distributed, do_rsgrid_replicated/),&
          default_i_val=do_rsgrid_replicated, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="EXTRAPOLATION",&
          variants=s2a("INTERPOLATION","WF_INTERPOLATION"),&
          description="Extrapolation strategy for e.g. MD",&
          usage="EXTRAPOLATION PS",&     ! XXXXXXXXXXX how to add the order ? i.e. PS 3
          enum_c_vals=s2a("PS","USE_GUESS","USE_PREV_P","USE_PREV_WF",&
          "USE_PREV_RHO_R", "LINEAR_WF", "LINEAR_P", "LINEAR_PS"),&
          enum_i_vals=(/ do_extrpl_ps,&
                         do_extrpl_use_guess,&
                         do_extrpl_use_prev_p,&
                         do_extrpl_use_prev_wf,&
                         do_extrpl_use_prev_rho_r,&
                         do_extrpl_linear_wf,&
                         do_extrpl_linear_p,&
                         do_extrpl_linear_ps /),&
          default_i_val=do_extrpl_linear_p, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="METHOD",&
          description="What method should be employed",&
          usage="METHOD GAPW",&
          enum_c_vals=s2a("GAPW","GPW","KG_POL","KG_NOPOL","KG_GPW"),&
          enum_desc=s2a("Gaussian and augmented plane waves method",&
                        "Gaussian and plane waves method",&
                        "Polarisable Kim-Gordon method",&
                        "Non-polarisable Kim-Gordon method", " ? "),&
          enum_i_vals=(/ do_method_gapw, do_method_gpw, do_method_kg_pol,&
                         do_method_kg_nopol, do_method_kg_gpw/),&
          default_i_val=do_method_gpw, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

  END IF

END SUBROUTINE create_qs_section

!!****f* input_cp2k/create_tddfpt_section [1.0] *
!!
!!   NAME
!!     create_tddfpt_section
!!
!!   FUNCTION
!!     creates the input section for the tddfpt part
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
!!     teo
!!
!!*** **********************************************************************
SUBROUTINE create_tddfpt_section(section,error)
  ! Arguments
  TYPE(section_type), pointer :: section
  type(cp_error_type), optional, intent(inout) :: error
  ! Local Variables
  logical :: failure
  character(len=*), parameter :: routineN='create_tddfpt_section',&
        routineP=moduleN//':'//routineN
  TYPE(section_type), pointer :: subsection
  TYPE(keyword_type), pointer :: keyword

  failure=.false.

  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL section_create(section,"tddfpt",&
          description="parameters needed to set up the Time Dependent Density Functional PT",&
          n_keywords=5, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
          error=error)

     NULLIFY(subsection,keyword)

     ! Integer
     CALL keyword_create(keyword, name="NEV",&
          variants=(/"N_EV",&
                     "EV  "/),&
          description=" number of excitations to calculate",&
          usage="NEV someInteger>0",&
          n_var=1,type_of_var=integer_t,&
          required=.TRUE., error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="NREORTHO",&
          variants=s2a("N_REORTHO","REORTHO","REORTHOGONALITAZIONS"),&
          description=" number of reorthogonalization steps",&
          usage="NREORTHO someInteger>0",&
          n_var=1,type_of_var=integer_t,&
          required=.TRUE., error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)    

     ! Logical
     CALL keyword_create(keyword, name="KERNEL",&
          description="compute the kernel (debug purpose only)",&
          usage="KERNEL",&
          default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     ! Character
     CALL keyword_create(keyword, name="RES_ETYPE",&
          variants=(/"RESTRICTED_EXCITATIONS_TYPE",&
                     "RES_E_TYPE                 "/),& 
          description="(singlets/triplets) for restricted calculation",&
          usage="RES_ETYPE T",&
          enum_c_vals=s2a("S","SINGLET","SINGLETS","T","TRIPLET","TRIPLETS"),&
          enum_i_vals=(/ do_tddfpt_singlet, do_tddfpt_singlet, do_tddfpt_singlet,&
                         do_tddfpt_triplet, do_tddfpt_triplet, do_tddfpt_triplet/),&
          default_i_val=do_tddfpt_singlet, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)    
     
     CALL keyword_create(keyword, name="DIAG_METHOD",&
          variants=(/"DIAGONALIZATION_METHOD"/),& 
          description="Diagonalization method used in tddfpt.",&
          usage="DIAG_METHOD DAVIDSON",&
          enum_c_vals=s2a("DAVIDSON","LANCZOS"),&
          enum_i_vals=(/ do_tddfpt_davidson, do_tddfpt_lanczos/),&
          default_i_val=do_tddfpt_lanczos, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)      

     CALL keyword_create(keyword, name="OE_CORR",&
          variants=(/"ORBITAL_EIGENVALUES_CORRECTION"/),& 
          description="Which type of orbital eigenvalue correction to use\n"//&
                      "(to yield better HOMO-LUMO energies)",&
          usage="OE_CORR SAOP",&
          enum_c_vals=s2a("NONE","SAOP","SIC"),&
          enum_i_vals=(/ oe_none, oe_saop, oe_sic /),&
          default_i_val=oe_none, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)      

     CALL create_xc_section(subsection,error)
     CALL section_add_subsection(section, subsection, error=error)
     CALL section_release(subsection,error=error)

  END IF


  END SUBROUTINE create_tddfpt_section

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

    CALL section_create(subsection, name="xc_grid",&!FM to do
         description="The xc parameters used when calculating the xc on the grid",&
         n_keywords=5, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
         error=error)
       CALL keyword_create(keyword, name="xc_smooth_rho",&
            description="The density smoothing used for the xc calculation",&
            usage="xc_smooth_rho nn10", default_i_val=xc_rho_no_smooth,&
            enum_c_vals=s2a("NONE","NN50","NN10","SPLINE2","NN6","SPLINE3","NN4"),&
            enum_i_vals=(/ xc_rho_no_smooth, xc_rho_nn50, xc_rho_nn10, &
            xc_rho_spline2_smooth,xc_rho_spline2_smooth,xc_rho_spline3_smooth,&
            xc_rho_spline3_smooth/),&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="xc_deriv",&
            description="The method used to compute the derivatives",&
            usage="xc_deriv NN10_SMOOTH", default_i_val=xc_deriv_spline2,&
            enum_c_vals=s2a("PW","SPLINE3","SPLINE2","NN50_SMOOTH","NN10_SMOOTH",&
            "SPLINE2_SMOOTH", "NN6_SMOOTH", "SPLINE3_SMOOTH", "NN4_SMOOTH"),&
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
    TYPE(section_type), pointer :: subsection

    failure=.false.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"scf",&
            description="parameters needed perform an scf run",&
            n_keywords=24, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       nullify(subsection)
       CALL create_ot_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

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
            enum_c_vals=s2a("ATOMIC","RESTART","RANDOM","CORE"),&
            enum_i_vals=(/atomic_guess,restart_guess,random_guess,core_guess/),&
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

    CALL create_qmmm_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

!FM    CALL create_metadyn_section(subsection,error)
!FM    CALL section_add_subsection(section, subsection, error=error)
!FM    CALL section_release(subsection,error=error)

    CALL create_poisson_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)    

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
            enum_c_vals=s2a("init","pos","posvel"),&
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
!!****f* input_cp2k/create_ewald_section *
!!
!!   NAME
!!     create_ewald_section
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: will contain the ewald section 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************

SUBROUTINE create_ewald_section(section,error)
  TYPE(section_type), pointer :: section
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='create_ewald_section',&
        routineP=moduleN//':'//routineN
  TYPE(keyword_type), pointer :: keyword

  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  if (.not. failure) then
    CALL section_create(section,name="ewald",&
         description="ewald parameters",&
         n_keywords=7, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
         error=error)

    nullify(keyword)
    CALL keyword_create(keyword, name="EWALD_TYPE",&
         description="The type of ewald you want to perform."//&
         " EWALD is the standard non-fft based ewald."//&
         " SPME is the smooth particle mesh using beta-Eurler splines."//&
         " PME is the particle mesh using fft interpolation.",&
         usage="EWALD_TYPE SPME}",&
         default_i_val=do_ewald,&
         enum_c_vals=s2a("EWALD", "SPME", "PME"),&
         enum_i_vals=(/do_ewald,do_spme,do_pme/),error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="UNITS",&
          variants=(/"UNIT"/),& 
          description="Specify the units used in the ewald input.",&
          usage="UNITS BOHR",&
          enum_c_vals=s2a("ANGSTROM","SCALED_ANGSTROM","BOHR","SCALED_BOHR"),&
          enum_i_vals=(/ use_unit_angstrom, use_unit_bohr,&
          use_unit_scaled_angstrom, use_unit_scaled_bohr/),&
          default_i_val=use_unit_angstrom, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)     


    CALL keyword_create(keyword, name="alpha",&
         description="alpha parameter associated with ewald",&
         usage="alpha .44", default_r_val=.44_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="gmax",&
         description="maximum number of grid points",&
         usage="gmax 21", default_i_val=21,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="ns_max",&
         description="number of grid points on small mesh (PME only)",&
         usage="ns_max 12", default_i_val=12,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="o_spline",&
         description="order of the beta-Euler spline",&
         usage="o_spline 6", default_i_val=6,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="epsilon",&
         description="tolerance of gaussians for fft interpolation",&
         usage="epsilon 1e-6", default_r_val=1.e-6_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

  end if
END SUBROUTINE create_ewald_section

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
    TYPE(keyword_type), pointer :: keyword
    TYPE(section_type), pointer :: subsection

    failure=.false.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       CALL section_create(section,name="subsys",&
            description="a subsystem: coordinates, topology, molecules and cell",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection,keyword)
       CALL create_cell_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="coord",&
            description="The coordinates for simple systems (like the QM ones)"//&
            " xyz format by default. More complex systems should be given with"//&
            " an external pdb file.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="The atoms coordinates in xyz format.",&
            usage="H 0.7 0.3 0.3", type_of_var=char_t, n_var=-2, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="kind",&
            description="The description of the kind of the atoms (mostly for QM)",&
            n_keywords=5, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)
       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="The name of the kind described in this section.",&
            usage="H", default_c_val="DEFAULT", error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AUX_BASIS_SET",&
            variants=s2a("AUXILIARY_BASIS_SET", "AUX_BASIS"),&
            description="The auxliliary basis set (GTO type)",&
            usage="AUX_BASIS_SET DZVP", type_of_var=char_t, &
            n_var=1, required=.TRUE., error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="AUX_BASIS_NORMALIZATION",&
            variants=s2a( "AUXILIARY_BASIS_NORMALIZATION", "AUX_BASIS_NORM"),&
            description="The normalization of the auxliliary basis set",&
            usage="AUX_BASIS_NORMALIZATION NO", &
            enum_c_vals=s2a( "NONE","NO","WFN","WAVEFUNCTION","FUNCTION",&
            "DENS","DENSITY","SQUARE","UNDEFINED"),&
            enum_i_vals=(/0,0,1,1,1,2,2,2,-1/),&
            default_i_val=-1,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="BASIS_SET",&
            variants=s2a("ORBITAL_BASIS_SET","ORB_BASIS"),&
            description="The basis set (GTO type)",&
            usage="BASIS_SET DZVP", type_of_var=char_t, &
            n_var=1, required=.TRUE., error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="BASIS_NORMALIZATION",&
            variants=s2a("ORBITAL_BASIS_NORMALIZATION","ORB_BASIS_NORM"),&
            description="The normalization of the auxliliary basis set",&
            usage="AUX_BASIS_NORMALIZATION NO", &
            enum_c_vals=s2a( "NONE", "NO", "WFN", "WAVEFUNCTION",&
            "FUNCTION", "DENS", "DENSITY", "SQUARE", "UNDEFINED"),&
            enum_i_vals=(/0,0,1,1,1,2,2,2,-1/),&
            default_i_val=-1,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="ELEMENT",&
            variants=(/ "ELEMENT_SYMBOL" /),&
            description="The element of the actual kind "//&
            "(if not given it is inferred from the kind name)",&
            usage="ELEMENT O", type_of_var=char_t,n_var=1,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="MASS",&
            variants=s2a("ATOMIC_MASS", "ATOMIC_WEIGHT", "WEIGHT"),&
            description="The mass of the atom "//&
            "(if negative or non present it is inferred from the element symbol)",&
            usage="MASS 2.0", type_of_var=real_t,n_var=1,&
            required=.FALSE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="POTENTIAL",&
            variants=(/ "POT" /),&
            description="The element of the actual kind "//&
            "(if not given it is inferred from the kind name)",&
            usage="POTENTIAL ghost", type_of_var=char_t, n_var=1,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="HARD_EXP_RADIUS",&
            description="The region where the hard density is supposed to be confined (GAPW)(in Bohr)",&
            usage="HARD_EXP_RADIUS 0.9", default_r_val=0.8_dp,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RHO0_EXP_RADIUS",&
            description="the radius which defines the atomic region where "//&
            "the hard compensation density is confined."//&
            "should be less than HARD_EXP_RADIUS (GAPW)[Bohr]",&
            usage="RHO_EXP_RADIUS 0.9", default_r_val=0.8_dp,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="LEBEDEV_GRID",&
            description="The number of points for the angular part of "//&
            "the local grid (GAPW)",&
            usage="LEBEDEV_GRID 40", default_i_val=50,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RADIAL_GRID",&
            description="The number of points for the radial part of "//&
            "the local grid (GAPW)",&
            usage="RADIAL_GRID 70", default_i_val=50,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

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
!!****f* input_cp2k/create_qmmm_section [1.0] *
!!
!!   NAME
!!     create_qmmm_section
!!
!!   FUNCTION
!!     Creates the QM/MM section
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
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_qmmm_section(section,error)
    TYPE(section_type), pointer :: section
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_qmmm_section',&
         routineP=moduleN//':'//routineN
    type(section_type), pointer :: subsection 
    TYPE(keyword_type), pointer :: keyword

    failure=.false.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="qmmm",&
            description="Input for QM/MM calculations.",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="A_QM",&
            description="Set the A vector of the QM Cell to the Reals value. The other values will be"//&
            "referred to the MM cell. It Changes the shape of the QM cell in one dimension.",&
            usage="A_QM real real real",&
            n_var=3,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="B_QM",&
            description="Set the B vector of the QM Cell to the Reals value. The other values will be"//&
            "referred to the MM cell. It Changes the shape of the QM cell in one dimension.",&
            usage="B_QM real real real",&
            n_var=3,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="C_QM",&
            description="Set the C vector of the QM Cell to the Reals value. The other values will be"//&
            "referred to the MM cell. It Changes the shape of the QM cell in one dimension.",&
            usage="C_QM real real real",&
            n_var=3,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ABC_QM",&
            description="Set the A,B,C vectors of  QM Cell to the Reals value."//&
            "It changes the overall shape of the QM cell.",&
            usage="ABC_QM real real real",&
            n_var=3,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="eps_mm_rspace",&
            description="Set the threshold g(r) - threshold = 0. Where g(r) is a gaussian function."//&
            "Standard value set to 1.0D-10",&
            usage="eps_mm_rspace real",&
            default_r_val=1.0E-10_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)     

       CALL keyword_create(keyword, name="E_COUPL",&
            variants=(/"QMMM_COUPLING",&
            "ECOUPL       "/),& 
            description="Set the Coupling's type between QM and MM system. Possible Options:"//&
            "(1) NONE: Just Mechanical Coupling. (2) GRID: Computing the electrostatic"//&
            "potential analytically on the Grid. (3) GAUSS: Mapping the MM electrostatic"//&
            "potential on the Grid with Gaussian.",&
            usage="E_COUPL char",&
            enum_c_vals=s2a( "NONE","GRID","GAUSS"),&
            enum_i_vals=(/ do_qmmm_none, do_qmmm_grid, do_qmmm_gauss/),&
            default_i_val=do_qmmm_none, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)    

       CALL section_create(subsection,name="QM_KIND",&
            description="Information about the kind qm that some atom have",&
            n_keywords=4, n_subsections=0, repeats=.true., required=.true.,&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="The qm kind",&
            usage="O",&
            n_var=1,type_of_var=char_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MM_INDEX",&
            description="The indexes of the mm atoms that have this kind",&
            usage="MM_INDEX 1 2",&
            n_var=-1,type_of_var=integer_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RADIUS",&
            description="Specifies the radius of the atomic kinds [angstrom]",&
            usage="RADIUS ",n_var=1,type_of_var=real_t,&
            default_r_val=0.8_dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SUBSYS",&
            description="Specifies which subsystem the qm atom belongs to.",&
            usage="SUBSYS ",n_var=1,type_of_var=integer_t,&
            default_i_val=1,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)    

    END IF

  END SUBROUTINE create_qmmm_section
!***************************************************************************
!!****f* input_cp2k/create_poisson_section [1.0] *
!!
!!   NAME
!!     create_poisson_section
!!
!!   FUNCTION
!!     Creates the Poisson section
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
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_poisson_section(section,error)
    TYPE(section_type), pointer :: section
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_poisson_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword

    failure=.false.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       CALL section_create(section,name="poisson",&
            description="Sets up the poisson resolutor.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       nullify(keyword)
       CALL keyword_create(keyword, name="POISSON_SOLVER",&
            variants=(/"POISSON", "PSOLVER"/),& 
            description="Specify which kind of solver to use to solve the Poisson equation.",&
            usage="POISSON_SOLVER char",&
            enum_c_vals=s2a( "PERIODIC3D", "ANALYTIC2D", "ANALYTIC1D", "ANALYTIC0D", "MT2D",&
            "MT1D", "MT0D", "HOCKNEY2D", "HOCKNEY1D", "HOCKNEY0D"),&
            enum_i_vals=(/ periodic3d, analytic2d, analytic1d, analytic0d, mt2d, &
            mt1d, mt0d, hockney2d, hockney1d, hockney0d/),&
            default_i_val=periodic3d, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)    

    end if
  END SUBROUTINE create_poisson_section
!***************************************************************************


!!****f* input_cp2k/create_cell_section *
!!
!!   NAME
!!     create_cell_section
!!
!!   FUNCTION
!!     creates the cell section
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_cell_section(section,error)
    TYPE(section_type), pointer :: section
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_cell_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword
    failure=.false.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"cell",&
            description="parameters needed to set up the cell",&
            n_keywords=6, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       nullify(keyword)
! Real
       CALL keyword_create(keyword, name="A",&
            description="Read the hmat(1,1); hmat(2,1); hmat(3,1)",&
            usage="A real real real",&
            n_var=-1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="B",&
            description="Read the hmat(1,2); hmat(2,2); hmat(3,2)",&
            usage="B real real real",&
            n_var=-1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="C",&
            description="Read the hmat(1,3); hmat(2,3); hmat(3,3)",&
            usage="C real real real",&
            n_var=-1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ABC",&
            description="Read the hmat(1,1); hmat(2,2); hmat(3,3)",&
            usage="ABC real real real",&
            n_var=-1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SUBCELLS",&
            description="??? SUBCELLS ???",&
            usage="SUBCELLS real",&
            n_var=-1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

! Character
       CALL keyword_create(keyword, name="UNITS",&
            variants=(/"UNIT"/),& 
            description="Specify the global units used in the input.",&
            usage="UNITS ANGSTROM",&
            enum_c_vals=s2a( "ANGSTROM","SCALED_ANGSTROM","BOHR","SCALED_BOHR"),&
            enum_i_vals=(/ use_unit_angstrom, use_unit_bohr,&
            use_unit_scaled_angstrom, use_unit_scaled_bohr/),&
            default_i_val=use_unit_bohr, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)     

    END IF
  END SUBROUTINE create_cell_section

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
       CALL keyword_create(keyword, name="REL_CUTOFF",&
            variants=(/"RELATIVE_CUTOFF"/),&
            description="Determines the grid at which a Gaussian is mapped,"//&
            " giving the cutoff used for a gaussian with alpha=1",&
            usage="RELATIVE_CUTOFF real", default_r_val=15.0_dp,&
            error=error)
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
            enum_c_vals=s2a("pw","linear","spline3_nopbc","spline3"),&
            enum_i_vals=(/pw_interp,linear_interp,&
            spline3_nopbc_interp,spline3_pbc_interp/),error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="aint_precond",&
            description="the approximate inverse to use to get the starting point"//&
            " for the linear solver of the spline3 methods",&
            usage="kind spline3",&
            default_i_val=pw_interp,&
            enum_c_vals=s2a( "copy","spl3_nopbc_aint","spl3_nopbc_precond1"),&
            enum_i_vals=(/no_precond,precond_spl3_aint, precond_spl3_1/),&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="precond",&
            description="The preconditioner used"//&
            " for the linear solver of the spline3 methods",&
            usage="kind spline3",&
            default_i_val=pw_interp,&
            enum_c_vals=s2a("copy","spl3_nopbc_aint","spl3_nopbc_precond1"),&
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

!!****f* input_cp2k/create_sic_section *
!!
!!   NAME
!!     create_sic_section
!!
!!   FUNCTION
!!     creates the sic (self interaction correction) section
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
SUBROUTINE create_sic_section(section,error)
  TYPE(section_type), pointer :: section
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='create_sic_section',&
        routineP=moduleN//':'//routineN
  TYPE(section_type), POINTER :: subsection,print_key
  TYPE(keyword_type), pointer :: keyword
  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL section_create(section,"sic",&
         description="parameters for the self interaction correction",&
         n_keywords=6, n_subsections=0, repeats=.FALSE., required=.false.,&
         error=error)

     nullify(keyword)

     CALL keyword_create(keyword, name="SIC_SCALING_A",&
          description="Scaling of the coulomb term in sic [experimental]",&
          usage="SIC_SCALING_A 0.5",&
          default_r_val=1.0_dp,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="SIC_SCALING_B",&
          description="Scaling of the xc term in sic [experimental]",&
          usage="SIC_SCALING_B 0.5",&
          default_r_val=1.0_dp,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="SIC_METHOD",&
          description="Method used to remove the self interaction",&
          usage="SIC_METHOD MAURI_US",&
          default_i_val=sic_none,&
          enum_c_vals=s2a( "NONE", "MAURI_US", "MAURI_SPZ", "AD"),&
          enum_i_vals=(/sic_none,sic_mauri_us,sic_mauri_spz,sic_ad/),&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

  END IF
END SUBROUTINE create_sic_section

!!****f* input_cp2k/create_ot_section *
!!
!!   NAME
!!     create_ot_section
!!   FUNCTION
!!     makes the orbital transformation section
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     11.2004 created [Joost VandeVondele]
!!
!!   SOURCE
!!*** **********************************************************************
SUBROUTINE create_ot_section(section,error)
  TYPE(section_type), pointer :: section
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='create_ot_section',&
        routineP=moduleN//':'//routineN
  TYPE(section_type), POINTER :: subsection,print_key
  TYPE(keyword_type), pointer :: keyword
  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL section_create(section,"ot",&
         description="Sets the various options for the orbital transformation (OT) method. "//&
                     "Default settings already provide an efficient, yet robust method. "//&
                     "Well-behaved systems might benefit from using a DIIS minimizer.",&
         n_keywords=6, n_subsections=0, repeats=.FALSE., required=.false.,&
         error=error)

     nullify(keyword)

     CALL keyword_create(keyword, name="MINIMIZER",&
          description="Minimizer to be used with the OT method",&
          usage="MINIMIZER DIIS",&
          default_i_val=ot_mini_cg,&
          enum_c_vals=s2a( "SD", "CG", "DIIS"),&
          enum_desc=s2a("Steepest descent","Conjugate Gradients", &
                        "Direct inversion in the iterative subspace"),&
          enum_i_vals=(/ot_mini_sd,ot_mini_cg,ot_mini_diis/),&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="SAFE_DIIS",&
            variants=(/"SAFER_DIIS"/),&
            description="Reject DIIS steps if they point away from the"//&
            "direction of minus, do SD in that case ",&
            usage="SAFE_DIIS ON", default_l_val=.true.,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="N_DIIS",&
          variants=(/"NDIIS"/),&
          description="Number of history vectors to be used with DIIS",&
          usage="N_DIIS 4",&
          default_i_val=7,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="LINESEARCH",&
          variants=(/"LINE_SEARCH"/),&
          description="1D line search algorithm to be used with the OT minimizer,"//&
                      " in increasing order of robustness and cost. MINIMIZER CG combined with"//&
                      " LINESEARCH GOLD should always find an electronic minimum.",&
          usage="LINESEARCH GOLD",&
          default_i_val=ot_ls_2pnt,&
          enum_c_vals=s2a( "NONE", "2PNT", "3PNT","GOLD"),&
          enum_desc=s2a("take fixed lenght steps","extrapolate based on 2 points", &
                        "... or on 3 points","perform 1D golden section search of the minimum"),&
          enum_i_vals=(/ot_ls_none,ot_ls_2pnt,ot_ls_3pnt,ot_ls_gold/),&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="STEPSIZE",&
          description="Initial stepsize used for the line search",&
          usage="STEPSIZE 0.4",&
          default_r_val=0.15_dp,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="GOLD_TARGET",&
          description="Target relative uncertainty in the minimum position of the minimum",&
          usage="GOLD_TARGET 0.1",&
          default_r_val=0.01_dp,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="PRECONDITIONER",&
          description="Type of preconditioner to be used with all minimization schemes."//&
                      "They differ in effectiveness, cost of construction, cost of application.",&
          usage="PRECONDITIONER FULL_SINGLE",&
          default_i_val=ot_precond_full_kinetic,&
          enum_c_vals=s2a("FULL_SINGLE","FULL_KINETIC","FULL_S_INVERSE","SPARSE_DIAG","NONE"),&
          enum_desc=s2a("Based on diagonalisation","Cholesky inversion of S and T", &
                        "Cholesky inversion of S","Only based on atomic blocks", &
                        "skip preconditioning"),&
          enum_i_vals=(/ot_precond_full_single,ot_precond_full_kinetic,ot_precond_s_inverse,&
                        ot_precond_sparse_diag,ot_precond_none/),&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="ENERGY_GAP",&
          description="Energy gap estimate [a.u.] for preconditioning",&
          usage="EPS_TAYLOR 0.1",&
          default_r_val=0.2_dp,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="EPS_TAYLOR",&
          variants=(/"EPSTAYLOR"/),&
          description="Target accuracy of the taylor expansion for the matrix functions",&
          usage="EPS_TAYLOR 1.0E-15",&
          default_r_val=1.0E-16_dp,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="MAX_TAYLOR",&
          description="Maximum order of the Taylor expansion before diagonalisation is prefered",&
          usage="MAX_TAYLOR 5",&
          default_i_val=4,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="ROTATION",&
          description="Introduce additional variables so that rotations of the occupied"//&
                      " subspace are allowed as well",&
          usage="ROTATION",lone_keyword_l_val=.TRUE.,&
          default_l_val=.false.,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

  END IF
END SUBROUTINE create_ot_section

END MODULE input_cp2k
