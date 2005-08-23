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
  USE kinds, ONLY: dp
  USE cp_para_types, ONLY: cp_para_env_type
  USE timings, ONLY: timeset, timestop
  USE input_section_types, ONLY: section_type, section_create, &
       section_add_keyword, section_add_subsection, section_release,&
       section_vals_type,&
       section_vals_get_subs_vals, section_vals_get
  USE input_constants
  USE input_keyword_types, ONLY: keyword_type, keyword_create,&
       keyword_release
  USE input_val_types, ONLY: val_type, no_t, logical_t, integer_t, real_t, char_t,&
       enum_t
  USE cp_output_handling, only: cp_print_key_section_create

  USE string_utilities, ONLY: s2a
  IMPLICIT NONE
  PRIVATE
  
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='input_cp2k'
  
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
       NULLIFY(section)

       CALL create_global_section(section,error=error)
       CALL section_add_subsection(root_section,section,error=error)
       CALL section_release(section,error=error)

       CALL create_motion_section(section,error=error)
       CALL section_add_subsection(root_section,section,error=error)
       CALL section_release(section,error=error)

       CALL create_force_eval_section(section,error=error)
       CALL section_add_subsection(root_section,section,error=error)
       CALL section_release(section,error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword,name="_DEFAULT_KEYWORD_",&
            description="ignores lines outside sections",&
            usage=" ", type_of_var=no_t,error=error)
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
  FUNCTION parsed_cp2k_input(input_file,check_this_section,error) RESULT(res)
    TYPE(section_vals_type), pointer :: input_file
    LOGICAL, INTENT(IN), OPTIONAL :: check_this_section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    logical :: res

    LOGICAL :: failure, my_check
    character(len=*), parameter :: routineN='parsed_cp2k_input',&
         routineP=moduleN//':'//routineN
    TYPE(section_vals_type), POINTER ::glob_section

    failure=.false.
    my_check = .FALSE.
    IF (PRESENT(check_this_section)) my_check = check_this_section
    res=ASSOCIATED(input_file)
    IF (res) THEN
       CPPrecondition(input_file%ref_count>0,cp_failure_level,routineP,error,failure)
       IF (.NOT.failure) THEN
          IF (.NOT.my_check) THEN
             glob_section => section_vals_get_subs_vals(input_file,"GLOBAL",&
                  error=error)
             CALL section_vals_get(glob_section,explicit=res,error=error)
          ELSE
             CALL section_vals_get(input_file,explicit=res,error=error)
          END IF
       END IF
    END IF
  END FUNCTION parsed_cp2k_input
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
         default_i_val=blacs_square,enum_c_vals=s2a("SQUARE","ROW","COLUMN"),&
         enum_desc=s2a("Distribution by marix blocks", "Distribution by matrix rows",&
         "Distribution by matrix columns"), &
         enum_i_vals=(/blacs_square,blacs_row,blacs_col/), error=error)
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
         enum_c_vals=program_names,&
         enum_desc=program_name_desc,&
         enum_i_vals=program_name_index,&
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

    CALL keyword_create(keyword, name="ERROR_FILE_NAME",&
         description="Name of the error file.",&
         usage="ERROR_FILE_NAME {filename}",&
         default_c_val="ERROR",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="RESTART_FILE_NAME",&
         description="Name of the restart file.",&
         usage="RESTART_FILE_NAME {filename}",&
         default_c_val="RESTART",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="RUN_TYPE",&
         description="Type of run that you want to perform Geometry "//&
         "optimization, md, montecarlo,...",&
         usage="RUN_TYPE MD",&
         default_i_val=energy_force_run,&
         enum_c_vals=run_type_names,&
         enum_i_vals=run_type_index,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="OUTPUT",&
         description="Specifies the format of the output file for geometry and velocity.",&
         usage="OUTPUT (XYZ|CRD|ATOMIC|XMOL)",&
         default_i_val=dump_xyz,&
         enum_c_vals=output_name_types,&
         enum_i_vals=output_name_index,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)    

     CALL keyword_create(keyword, name="WALLTIME",&
          variants=(/"WALLTI"/),&
          description="Maximum execution time for this run. Time in seconds.",&
          usage="WALLTIME real",  default_r_val=-1.0_dp, error=error)
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
    CALL keyword_create(keyword, name="PP_LIBRARY_PATH",&
         description="Path of the directory where PseudoPotential are located.",&
         usage="PP_LIBRARY_PATH path",&
         default_c_val="./",error=error)
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

    CALL create_poisson_dft_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_xas_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_print_dft_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

  END SUBROUTINE create_dft_section
!***************************************************************************


!****f* input_cp2k/create_localize_section *
!!
!!   NAME
!!     create_localize_section
!!   FUNCTION
!!     parameters fo the localization of wavefunctions
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     03.2005 created [MI]
!!
!!   SOURCE
!!*** **********************************************************************
  SUBROUTINE create_localize_section(section, error)

    TYPE(section_type), pointer :: section
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_localize_section',&
          routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword
    failure=.false.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

      NULLIFY(keyword)
      CALL section_create(section,name="LOCALIZE",&
           description="Use one of the available methods to define the localization "//&
           " and possibly to optimize it to a minimum or a maximum.",&
           n_keywords=8, n_subsections=0, repeats=.TRUE., required=.false.,&
           error=error)

      CALL keyword_create(keyword,name="MIN_OR_MAX",&
           description="Requires the maximization of the spread of the wfn",&
           usage="MIN_OR_MAX (MINIMIZE|MAXIMIZE)",&
           enum_c_vals=(/"SPREADMIN","SPREADMAX"/),&
           enum_i_vals=(/do_loc_min, do_loc_max/),&
           default_i_val=do_loc_min,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)
 
 
      CALL keyword_create(keyword,name="METHOD",&
           description="Method of optimization if any",&
           usage="METHOD (JACOBI|CRAZY|DIRECT|NONE)",&
           enum_c_vals=s2a("NONE","JACOBI","CRAZY","DIRECT"),&
           enum_i_vals=(/do_loc_jacobi, do_loc_crazy, do_loc_direct, do_loc_none/),&
           default_i_val=do_loc_jacobi,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)
 
      CALL keyword_create(keyword,name="OPERATOR",&
           description="Type of opertator which defines the spread functional",&
           usage="OPERATOR (BERRY|BOYS|PIPEK)",&
           enum_c_vals=s2a("BERRY","BOYS","PIPEK"),&
           enum_i_vals=(/op_loc_berry, op_loc_boys, op_loc_pipek/),&
           default_i_val=op_loc_berry,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)
 
      CALL keyword_create(keyword,name="SET_OF_STATES",&
           description="Selected states to be localized",&
           usage="SET_OF_STATES (ALL|UPTON|FROMITOJ|LIST|NONE)",&
           enum_c_vals=s2a("ALL","UPTON","FROMITOJ","LIST","NONE"),&
           enum_i_vals=(/state_loc_all, state_loc_upton, state_loc_fromitoj,& 
           state_loc_list,state_loc_none/),&
           default_i_val=state_loc_all,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)


      CALL keyword_create(keyword,name="LU_BOUNDS",&
           description="The lower and upper index of the set of wfn to be localized",&
           usage="LU_BOUNDS integer integer",&
           n_var=2,default_i_vals=(/0,0/),type_of_var=integer_t,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,name="UP_BOUND",&
           description="Highest index of  wfn to be localized",&
           usage="UP_BOUND  integer",&
           default_i_val=0,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)


      CALL keyword_create(keyword,name="LISTSTATES",&
           description="Indexes of the wfn to be localized"//&
           "This keyword can be repeated several times"//&
           "(useful if you have to specify many indexes).",&
           usage="LISTSTATES 1 2",&
           n_var=-1,type_of_var=integer_t,repeats=.TRUE.,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)
 
      CALL keyword_create(keyword,name="EACH_STEP",&
           description="Evry how many steps the localization calculation"//&
           " is performed  along an iterative procedure (e.g. MD)",&
           usage="EACH_STEP integer",&
           default_i_val=1,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="FNAME",&
           description="Root of the file names where to store the spectra",&
           usage="FNAME loc_orb",&
           default_c_val="loc_orb",error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_localize_section

!***************************************************************************
!!****f* input_cp2k/create_print_dft_section [1.0] *
!!
!!   NAME
!!     create_print_dft_section
!!
!!   FUNCTION
!!     Create the print dft section
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
  SUBROUTINE create_print_dft_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_print_dft_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: print_key
    TYPE(keyword_type), POINTER :: keyword
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="print",&
            description="Section of possible print options in DFT code.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(print_key, keyword)

       CALL cp_print_key_section_create(print_key,"ATOMIC_COORDINATES",&
            description="Controls the printing of atomic coordinates in the DFT code.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"BASIC_DATA_TYPES",&
            description="Controls the printing of the basic data types.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"BASIC_SETS",&
            description="Controls the printing of the basis sets used in the DFT calculation.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"BLACS_INFORMATION",&
            description="Controls the printing of information about the specified BLACS context.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"CELL_PARAMETERS",&
            description="Controls the printing of cell information", &
            print_level=silent_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"CORE_HAMILTONIAN_MATRIX",&
            description="Controls the printing of the core hamiltonian matrix.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"CORE_CHARGE_RADII",&
            description="Controls the printing of the core charge radii.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"CPUTIME",&
            description="Controls the printing of cputime information along the MM execution.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"DENSITY_MATRIX",&
            description="Controls the printing of the Density Matrix", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"DERIVATIVES",&
            description="Controls the printing of derivatives.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"DIIS_INFORMATION",&
            description="Controls the printing of the DIIS information.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"DISTRIBUTION",&
            description="Controls the printing of the distribution of the simple pair neighbour list.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"FORCE",&
            description="Controls the printing of the forces.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"HARTREE_MATRIX",&
            description="Controls the printing of the Hartree Matrix.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"INTERATOMIC_DISTANCES",&
            description="Controls the printing of interatomic distances.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"KIND_INFORMATION",&
            description="Controls the printing of information on atomic kinds.", &
            print_level=low_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"KIND_RADII",&
            description="Controls the printing of the kind radii", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"KINETIC_ENERGY_MATRIX",&
            description="Controls the printing the kinetic energy matrix.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"KS_MATRIX",&
            description="Controls the printing of the Kohn-Sham Matrix.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"MEMORY",&
            description="Controls the printing of information on memory usage", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"MO_EIGENVALUES",&
            description="Controls the printing of the eigenvalues of the MOs", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"MO_EIGENVECTORS",&
            description="Controls the printing of the Molecular Orbitals", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"MO_OCCUPATION_NUMBERS",&
            description="Controls the printing of the occupation numbers of MOs.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"NEIGHBOUR_LISTS",&
            description="Controls the printing of the set of neighbor lists.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"OPTICAL_CONDUCTIVITY",&
            description="Controls the printing of the optical conductivity.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"OVERLAP_MATRIX",&
            description="Controls the printing of the Overlap Matrix", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PGF_RADII",&
            description="Controls the printing of the radii of the primitive gaussian functions.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PHYSICAL_CONSTANTS",&
            description="Controls the printing of physical constants used in the code.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"POTENTIAL",&
            description="Controls the printing of the potential generated by the QM "//&
            "charge distribution on a cube file.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"POTENTIAL_ENERGY_MATRIX",&
            description="Controls the printing of the potential energy matrix", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PPL_RADII",&
            description="Controls the printing of the radii of the exponential functions "//&
            "of the Goedecker pseudopotential (GTH, local part)", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PPNL_RADII",&
            description="Controls the printing of the radii of the projector functions "//&
            "of the Goedecker pseudopotential (GTH, non-local part)", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PRJ_RADII",&
            description="Controls the printing of the radii of the one center projector", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PRECISION",&
            description="Controls the precision of the results printed.",&
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PROGRAM_BANNER",&
            description="Controls the printing of the banner of the MM program",&
            print_level=silent_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="Controls the printing of information regarding the run.",&
            print_level=silent_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PW_GRID_INFORMATION",&
            description="Controls the printing of information regarding the pw_grid structure.",&
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"RADII",&
            description="Controls the printing of the full species of radii", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SAB_ORB_NEIGHBOR_LIST",&
            description="Controls the printing of the orbital-orbital overlap neighbor lists", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SAC_PPL_NEIGHBOR_LIST",&
            description="Controls the printing of the orbital GTH-PPL operator overlap list.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SAC_PPNL_NEIGHBOR_LIST",&
            description="Controls the printing of the orbital GTH-PPNL operator overlap list.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SBC_PPNL_NEIGHBOR_LIST",&
            description="Controls the printing of the orbital GTH-PPNL operator overlap list.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SAC_OCE_NEIGHBOR_LIST",&
            description="Controls the printing of then orbital orbital-gapw projector  overlap list", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SBC_OCE_NEIGHBOR_LIST",&
            description="Controls the printing of then orbital orbital-gapw projector  overlap list", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SCF",&
            description="Controls the printing of the SCF parameters.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SET_RADII",&
            description="Controls the printing of the shell set radii", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SPHERICAL_HARMONICS",&
            description="Controls the printing of the Spherical Harmonics and the"//&
            "orbital transformation matrices", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"TOTAL_NUMBERS",&
            description="Controls the printing of the total number of kinds, atoms, etc.",&
            print_level=low_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)


       CALL cp_print_key_section_create(print_key,"W-MATRIX",&
            description="Controls the printing of the W matrix, involved in the force calculation.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"ORBITALS",&
            description="Controls the printing of the MOs with cube format", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)

       CALL keyword_create(keyword, name="NUMBER_ORBITALS",&
            description="Specify the number of orbitals to print (counting from HOMO to the lowest orbital)",&
            usage="NUMBER_ORBITALS {integer}",&
            default_i_val=5,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)       

       CALL keyword_create(keyword, name="ORBITAL_NUMBER",&
            description="Specify the orbital number to print.",&
            usage="ORBITAL_NUMBER {integer} .. {integer}",&
            required=.TRUE., type_of_var=integer_t, n_var=-1, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)       

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"E_DENSITY_CUBE",&
            description="Controls the printing of the electron density with cube format", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"V_HARTREE_CUBE",&
            description="Controls the printing of the hartree potential with a cube format", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"LUMO",&
            description="Controls the printing of the Lumo with cube format", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"HOMO",&
            description="Controls the printing of the Homo with cube format", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"LOCALIZATION",&
            Description="Controls the printing of the Wannier centers.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)

       CALL keyword_create(keyword, name="LOC_TYPE",&
            description="Specify the method to use for the localization.",&
            usage="LOC_TYPE (JACOBI|CRAZY|DIRECT|NONE)",&
            enum_c_vals=s2a("JACOBI","CRAZY","DIRECT","NONE"),&
            enum_i_vals=(/ do_loc_jacobi, do_loc_crazy, do_loc_direct, do_loc_none/),&
            default_i_val=do_loc_jacobi, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"MOLECULAR_STATES",&
            description="Controls the printing of the molecular states", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"HPM",&
            description="Controls the printing of the hpm info (if linked as external library)", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"DIPOLE",&
            description="Calculates the dipole, for the moment it only supports"//&
            "the stupid non periodic dipole, should be extended soon", &
            print_level=high_print_level,each=(/1/),filename=" ",&
            common_iter_levels=1,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"MULLIKEN",&
            description="Controls the printing of the Mulliken Population analysis.", &
            print_level=silent_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"LOCALIZED_WFN",&
            description="Controls the printing of new localized wavefunctions method.", &
            print_level=silent_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_print_dft_section


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
  TYPE(section_type), pointer :: subsection, print_key

  failure=.false.

  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL section_create(section,"qs",&
          description="parameters needed to set up the Quickstep framework",&
          n_keywords=34, n_subsections=2, repeats=.FALSE., required=.FALSE.,&
          error=error)

     NULLIFY(keyword, subsection)
    
    CALL keyword_create(keyword, name="RESTART_FILE_NAME",&
         description="The path to the pseudo potential file",&
         usage="RESTART_FILE_NAME ../RESTART",&
         default_c_val="RESTART",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
     ! Reals
     CALL keyword_create(keyword, name="EPS_DEFAULT",&
          description="Try setting all EPS_xxx to values leading to an energy correct up to EPS_DEFAULT",&
          usage="EPS_DEFAULT real", type_of_var=real_t,&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)     

     CALL keyword_create(keyword, name="EPS_CORE_CHARGE",&
          description="Precision for mapping the core charges",&
          usage="EPS_CORE_CHARGE real", default_r_val=1.0E-12_dp,&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)     

     CALL keyword_create(keyword, name="EPS_GVG_RSPACE",&
          variants=(/"EPS_GVG"/),&
          description="Sets precision of the realspace KS matrix element integration",&
          usage="EPS_GVG_RSPACE real",default_r_val=1.0E-5_dp ,&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="EPS_PGF_ORB",&
          description="Sets precision of the overlap matrix elements",&
          usage="EPS_PGF_ORB real",default_r_val=  1.0E-5_dp,&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)


     CALL keyword_create(keyword, name="EPS_PPL",&
          description="Sets precision of the local part of the pseudo potential",&
          usage="EPS_PPL real", default_r_val=1.0E-12_dp,&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="EPS_PPNL",&
          description="Sets precision of the non-local part of the pseudo potential",&
          usage="EPS_PPNL real",  default_r_val=1.0E-12_dp,&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="EPS_RHO",&
          description="Sets precision of the density mapping on the grids",&
          usage="EPS_RHO real",type_of_var=real_t,&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="EPS_RHO_RSPACE",&
          description="Sets precision of the density mapping in rspace",&
          usage="EPS_RHO_RSPACE real",default_r_val=1.0E-10_dp,&
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

     CALL keyword_create(keyword, name="LOCALIZE_WFN",&
          description="Check and possibly optimize the localization of the wavefunctions",&
          usage="LOCALIZE_WFN",&
          default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error) 

     NULLIFY(subsection)
     CALL create_localize_section(subsection,error)
     CALL section_add_subsection(section, subsection, error=error)
     CALL section_release(subsection,error=error)


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
          enum_c_vals=s2a("USE_GUESS","USE_PREV_P","USE_PREV_RHO_R","LINEAR_WF","LINEAR_P","LINEAR_PS","USE_PREV_WF","PS"),&
          enum_i_vals=(/ do_extrpl_use_guess,& 
                         do_extrpl_use_prev_p,&
                         do_extrpl_use_prev_rho_r,&
                         do_extrpl_linear_wf,&
                         do_extrpl_linear_p,&
                         do_extrpl_linear_ps,&
                         do_extrpl_use_prev_wf,&
                         do_extrpl_ps /),&
          default_i_val=do_extrpl_linear_p, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="EXTRAPOLATION_ORDER",&
          description="Order for the PS extrapolation",&
          usage="EXTRAPOLATION_ORDER {integer}",default_i_val=0, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)



     CALL keyword_create(keyword, name="METHOD",&
          description="What method should be employed",&
          usage="METHOD GAPW",&
          enum_c_vals=s2a("GAPW","GPW","KG_POL","KG_NOPOL","KG_GPW","EHT","MNDO","AM1","PM3"),&
          enum_desc=s2a("Gaussian and augmented plane waves method",&
                        "Gaussian and plane waves method",&
                        "Polarisable Kim-Gordon method",&
                        "Non-polarisable Kim-Gordon method", " DFT Embedding method ",&
                        "EHT semiempirical","MNDO semiempirical","AM1 semiempirical",&
                        "PM3 semiempirical"),&
          enum_i_vals=(/ do_method_gapw, do_method_gpw, do_method_kg_pol,&
                         do_method_kg_nopol, do_method_kg_gpw, do_method_eht,&
                         do_method_mndo, do_method_am1, do_method_pm3/),&
          default_i_val=do_method_gpw, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     NULLIFY(subsection)
     CALL create_se_control_section(subsection,error)
     CALL section_add_subsection(section, subsection, error=error)
     CALL section_release(subsection,error=error)

     NULLIFY(subsection)
     CALL create_mulliken_section(subsection,error)
     CALL section_add_subsection(section, subsection, error=error)
     CALL section_release(subsection,error=error)

     NULLIFY(subsection)
     CALL create_harris_section(subsection,error)
     CALL section_add_subsection(section, subsection, error=error)
     CALL section_release(subsection,error=error) 

  END IF

END SUBROUTINE create_qs_section

SUBROUTINE create_se_control_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_se_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    NULLIFY(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      CALL section_create(section,"SE",&
           description="parameters needed to set up the Semiempirical method",&
           n_keywords=7, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
           error=error)

      CALL keyword_create(keyword, name="orthogonal_basis",&
           description="Use orthogonal basis",&
           usage="ORTHOGONAL_BASIS ",default_l_val=.FALSE., error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="EWALD",&
          description="Do ewald",&
           usage="EWALD ",default_l_val=.FALSE., error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="DELTA",&
          description="??",&
          usage="DELTA {real} ",default_r_val=1.e-6_dp, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="ALPHA_GAUSS",&
          description="??",&
          usage="ALPHA_GAUSS {real} ",default_r_val=1._dp, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="EPS_GAUSS",&
          description="??",&
          usage="EPS_GAUSS {real} ",default_r_val=1.e-8_dp, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="RC_EXCHANGE",&
          description="??",&
          usage="RC_EXCHANGE {real} ",default_r_val=25._dp, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)
 
     CALL keyword_create(keyword, name="RC_COULOMB",&
          description="??",&
          usage="RC_COULOMB {real} ",default_r_val=10._dp, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

    END IF

END SUBROUTINE create_se_control_section

SUBROUTINE create_mulliken_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mulliken_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    nullify(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      CALL section_create(section,"MULLIKEN",&
           description="parameters needed to calculate the Mulliken charges",&
           n_keywords=7, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
           error=error)

      CALL keyword_create(keyword, name="STRENGTH",&
           description="??",&
           usage="STRENGTH {real} ",default_r_val=0.1_dp, error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)
  
      CALL keyword_create(keyword, name="TARGET",&
           description="??",&
           usage="TARGET {real} ",default_r_val=1._dp, error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="NATOMS",&
           description="??",&
           usage="NATOMS {integer} ",default_i_val=0, error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="ATOMS",&
           description="Specifies a list of atoms to ....",&
           usage="LIST {integer} {integer} .. {integer}", required=.TRUE.,&
           n_var=-1, type_of_var=integer_t, error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)
    END IF

END SUBROUTINE create_mulliken_section


!!****f* xc_derivatives/create_xc_fun_section [1.0] *
!!   NAME
!!     create_xc_fun_section
!!
!!   FUNCTION
!!     creates the structure of the section needed to select the xc functional
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section that will be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE create_xc_fun_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_xc_fun_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

  failure=.FALSE.
  
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL section_create(section,name="xc_functional",&
         description="The xc functional to use",&
         n_keywords=0, n_subsections=4, repeats=.FALSE., required=.FALSE.,&
         error=error)
     
     NULLIFY(subsection,keyword)
     CALL section_create(subsection,name="BECKE88",&
          description="Uses the Becke 88 exchange functional",&
          n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
          error=error)
     CALL section_add_subsection(section,subsection,error=error)
     CALL section_release(subsection,error=error)

     CALL section_create(subsection,name="LYP",&
          description="Uses the LYP functional",&
          n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
          error=error)
     CALL section_add_subsection(section,subsection,error=error)
     CALL section_release(subsection,error=error)

     CALL section_create(subsection,name="PADE",&
          description="Uses the PADE functional",&
          n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
          error=error)
     CALL section_add_subsection(section,subsection,error=error)
     CALL section_release(subsection,error=error)

     CALL section_create(subsection,name="HCTH",&
          description="Uses the HCTH class of functionals",&
          n_keywords=0, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
          error=error)
     CALL keyword_create(keyword,name="PARAMETER_SET",&
          description="Which version of the parameters should be used",&
          usage="PARAMETER_SET 407",&
          enum_c_vals=(/"93 ","120","147","407"/),&
          enum_i_vals=(/93,120,147,407/),&
          default_i_val=120,error=error)
     CALL section_add_keyword(subsection,keyword,error=error)
     CALL keyword_release(keyword,error=error)
     CALL section_add_subsection(section,subsection,error=error)
     CALL section_release(subsection,error=error)

     CALL section_create(subsection,name="OPTX",&
          description="Uses the OPTX functional",&
          n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
          error=error)
     CALL section_add_subsection(section,subsection,error=error)
     CALL section_release(subsection,error=error)

     CALL section_create(subsection,name="CS1",&
          description="Uses the CS1 functional",&
          n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
          error=error)
     CALL section_add_subsection(section,subsection,error=error)
     CALL section_release(subsection,error=error)

     CALL section_create(subsection,name="XGGA",&
          description="Uses one of the XGGA functionals (optimized versions of "//&
          "some of these functionals might be available outside this section).",&
          n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
          error=error)
     CALL keyword_create(keyword,name="FUNCTIONAL",&
          description="Which one of the XGGA functionals should be used",&
          usage="FUNCTIONAL PW86X",&
          enum_c_vals=(/&
          "BECKE88X",&
          "PW86X   ",&
          "PW91X   ",&
          "PBEX    ",&
          "REV_PBEX",&
          "OPTX    "/),&
          enum_i_vals=(/xgga_b88x,xgga_pw86,xgga_pw91,xgga_pbex,xgga_revpbe,xgga_opt/),&
          default_i_val=xgga_b88x,error=error)
     CALL section_add_keyword(subsection,keyword,error=error)
     CALL keyword_release(keyword,error=error)
     CALL section_add_subsection(section,subsection,error=error)
     CALL section_release(subsection,error=error)

     CALL section_create(subsection,name="KE_GGA",&
          description="Uses one of the KE_GGA functionals (optimized versions of "//&
          "some of these functionals might be available outside this section).",&
          n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
          error=error)
     CALL keyword_create(keyword,name="FUNCTIONAL",&
          description="Which one of the KE_GGA functionals should be used",&
          usage="FUNCTIONAL LLP",&
          enum_c_vals=(/"OL1 ","OL2 ","LLP ","PW86","PW91","LC  ","T92 ","PBE "/),&
          enum_i_vals=(/ke_ol1,ke_ol2,ke_llp,ke_pw86,ke_pw91,ke_lc,ke_t92,ke_pbe/),&
          default_i_val=ke_ol1,error=error)
     CALL section_add_keyword(subsection,keyword,error=error)
     CALL keyword_release(keyword,error=error)
     CALL section_add_subsection(section,subsection,error=error)
     CALL section_release(subsection,error=error)

     CALL section_create(subsection,name="P86C",&
          description="Uses the P86C functional",&
          n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
          error=error)
     CALL section_add_subsection(section,subsection,error=error)
     CALL section_release(subsection,error=error)

     CALL section_create(subsection,name="PW92",&
          description="Uses the PerdewWang correlation functional.",&
          n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
          error=error)
     CALL keyword_create(keyword,name="PARAMETRIZATION",&
          description="Which one of parametrizations should be used",&
          usage="PARAMETRIZATION DMC",&
          enum_c_vals=(/&
          "ORIGINAL",&
          "DMC     ",&
          "VMC     "/),&
          enum_i_vals=(/c_pw92,c_pw92dmc,c_pw92vmc/),&
          default_i_val=c_pw92,error=error)
     CALL section_add_keyword(subsection,keyword,error=error)
     CALL keyword_release(keyword,error=error)
     CALL section_add_subsection(section,subsection,error=error)
     CALL section_release(subsection,error=error)

     CALL section_create(subsection,name="PZ81",&
          description="Uses the PZ functional.",&
          n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
          error=error)
     CALL keyword_create(keyword,name="PARAMETRIZATION",&
          description="Which one of parametrizations should be used",&
          usage="PARAMETRIZATION DMC",&
          enum_c_vals=(/&
          "ORIGINAL",&
          "DMC     ",&
          "VMC     "/),&
          enum_i_vals=(/c_pz,c_pzdmc,c_pzvmc/),&
          default_i_val=pz_orig,error=error)
     CALL section_add_keyword(subsection,keyword,error=error)
     CALL keyword_release(keyword,error=error)
     CALL section_add_subsection(section,subsection,error=error)
     CALL section_release(subsection,error=error)

     CALL section_create(subsection,name="TFW",&
          description="Uses the TFW functional",&
          n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
          error=error)
     CALL section_add_subsection(section,subsection,error=error)
     CALL section_release(subsection,error=error)
     
     CALL section_create(subsection,name="TF",&
          description="Uses the TF functional",&
          n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
          error=error)
     CALL section_add_subsection(section,subsection,error=error)
     CALL section_release(subsection,error=error)
     
     CALL section_create(subsection,name="VWN",&
          description="Uses the VWN functional",&
          n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
          error=error)
     CALL section_add_subsection(section,subsection,error=error)
     CALL section_release(subsection,error=error)
     
     CALL section_create(subsection,name="XALPHA",&
          description="Uses the XALPHA (SLATER) functional.",&
          n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
          error=error)
     CALL keyword_create(keyword,name="XA",&
          description="Value of the xa parameter (this does not change the exponent, "//&
          "just the mixing)",&
          usage="XA 0.7", default_r_val=2._dp/3._dp,error=error)
     CALL section_add_keyword(subsection,keyword,error=error)
     CALL keyword_release(keyword,error=error)
     CALL section_add_subsection(section,subsection,error=error)
     CALL section_release(subsection,error=error)
     
     CALL section_create(subsection,name="TPSS",&
          description="Uses the TPSS functional",&
          n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
          error=error)
     CALL section_add_subsection(section,subsection,error=error)
     CALL section_release(subsection,error=error)
     
     CALL section_create(subsection,name="PBEX",&
          description="Uses the PBE functional",&
          n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
          error=error)
     CALL section_add_subsection(section,subsection,error=error)
     CALL section_release(subsection,error=error)
          
     CALL section_create(subsection,name="PBEC",&
          description="Uses the PBE functional",&
          n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
          error=error)
     CALL section_add_subsection(section,subsection,error=error)
     CALL section_release(subsection,error=error)
          
  END IF
END SUBROUTINE create_xc_fun_section

!***************************************************************************
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
     CALL keyword_create(keyword, name="MAX_KV",&
          variants=s2a("MAX_VECTORS"),&
          description=" maximal number of Krylov space vectors",&
          usage="MAX_KV someInteger>0",&
          n_var=1,type_of_var=integer_t,&
          default_i_val=60, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="RESTARTS",&
          variants=s2a("N_RESTARTS"),&
          description=" maximal number subspace search restarts",&
          usage="RESTARTS someInteger>0",&
          n_var=1,type_of_var=integer_t,&
          default_i_val=5, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="NEV",&
          variants=s2a("N_EV", "EV"),&
          description=" number of excitations to calculate",&
          usage="NEV someInteger>0",&
          n_var=1,type_of_var=integer_t,&
          default_i_val=1, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="NREORTHO",&
          variants=s2a("N_REORTHO","REORTHO","REORTHOGONALITAZIONS"),&
          description=" number of reorthogonalization steps",&
          usage="NREORTHO someInteger>0",&
          n_var=1,type_of_var=integer_t,&
          default_i_val=2, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)    

     ! Logical
     CALL keyword_create(keyword, name="KERNEL",&
          variants=s2a("DO_KERNEL"),&
          description="compute the kernel (debug purpose only)",&
          usage="KERNEL logical_value",&
          default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="LSD_SINGLETS",&
          description="compute singlets using lsd vxc kernel",&
          usage="LSD_SINGLETS logical_value",&
          default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="INVERT_S",&
          variants=s2a("INVERT_OVERLAP"),&
          description="use the inverse of the overlap matrix",&
          usage="INVERT_S logical_value",&
          default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="PRECONDITIONER",&
          variants=s2a("PRECOND"),&
          description="use the preconditioner (only for Davidson)",&
          usage="PRECONDITIONER logical_value",&
          default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     ! Character
     CALL keyword_create(keyword, name="RES_ETYPE",&
          variants=s2a("RESTRICTED_EXCITATIONS_TYPE", "RES_E_TYPE"),& 
          description="(singlets/triplets) for restricted calculation",&
          usage="RES_ETYPE T",&
          enum_c_vals=s2a("S","SINGLET","SINGLETS","T","TRIPLET","TRIPLETS"),&
          enum_i_vals=(/ tddfpt_singlet, tddfpt_singlet, tddfpt_singlet,&
                         tddfpt_triplet, tddfpt_triplet, tddfpt_triplet/),&
          default_i_val=tddfpt_singlet, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)    
     
     CALL keyword_create(keyword, name="DIAG_METHOD",&
          variants=s2a("DIAGONALIZATION_METHOD", "METHOD"),& 
          description="Diagonalization method used in tddfpt",&
          usage="DIAG_METHOD DAVIDSON",&
          enum_c_vals=s2a("DAVIDSON","LANCZOS"),&
          enum_i_vals=(/ tddfpt_davidson, tddfpt_lanczos/),&
          default_i_val=tddfpt_davidson, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)      

     CALL keyword_create(keyword, name="OE_CORR",&
          variants=s2a("ORBITAL_EIGENVALUES_CORRECTION"),& 
          description="Which type of orbital eigenvalue correction to use\n"//&
                      "(to yield better HOMO-LUMO energies)",&
          usage="OE_CORR SAOP",&
          enum_c_vals=s2a("NONE", "LB", "LB_ALPHA", "LB94", "GLLB", "GLB", "SAOP","SIC"),&
          enum_i_vals=(/ oe_none, oe_lb, oe_lb, oe_lb, oe_gllb, oe_gllb, oe_saop, oe_sic /),&
          default_i_val=oe_none, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)      

     ! Real
     CALL keyword_create(keyword, name="CONVERGENCE",&
          variants=s2a("CONV"),&
          description="The convergence of the eigenvalues",&
          usage="CONVERGENCE 1.0E-6 ",&
          n_var=1,type_of_var=real_t,&
          default_r_val=1.0e-5_dp, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL create_xc_section(subsection,error)
     CALL section_add_subsection(section, subsection, error=error)
     CALL section_release(subsection,error=error)

     CALL create_sic_section(subsection,error)
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

       CALL keyword_create(keyword, name="DENSITY_SMOOTH_CUTOFF_RANGE",&
            description="Parameter for the smoothing procedure in"//&
            "xc calculation",&
            usage="gradient_cutoff {real}", default_r_val=0.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword, name="tau_cutoff",&
            description="The cutoff on tau used by the xc calculation",&
            usage="tau_cutoff 1.e-11", default_r_val=1.0e-10_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FUNCTIONAL_ROUTINE",&
       description="Select the code for xc calculation",&
            usage="FUNCTIONAL_ROUTINE NEW", default_i_val=xc_new_f_routine,&
            enum_c_vals=s2a("NEW","OLD","TEST_LSD","DEBUG"),&
            enum_i_vals=(/ xc_new_f_routine, xc_old_f_routine, &
            xc_test_lsd_f_routine, xc_debug_new_routine/),error=error)
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
            usage="xc_deriv NN10_SMOOTH", default_i_val=xc_deriv_pw,&
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
            usage="ADDED_MOS", default_i_val=0,n_var=-1,error=error)
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

    CALL create_mm_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_qmmm_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_bsse_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_subsys_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_properties_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

  END SUBROUTINE create_force_eval_section
!***************************************************************************


!***************************************************************************
!!****f* input_cp2k/create_properties_section [1.0] *
!!
!!   NAME
!!     create_properties_section
!!
!!   FUNCTION
!!     Create the PROPERTIES section 
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
  SUBROUTINE create_properties_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_properties_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    TYPE(section_type), POINTER :: subsection
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="PROPERTIES",&
            description="This section is used to set up the PROPERTIES calculation.",&
            n_keywords=3, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword, subsection)
    
       CALL create_linres_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_fitcharge_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
      
    END IF

  END SUBROUTINE create_properties_section

!***************************************************************************
!!****f* input_cp2k/create_fitcharge_section [1.0] *
!!
!!   NAME
!!     create_fitcharge_section
!!
!!   FUNCTION
!!     Create the BSSE section for counterpoise correction
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
  SUBROUTINE create_fitcharge_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_fitcharge_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    TYPE(section_type), POINTER :: subsection
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="FIT_CHARGE",&
            description="This section is used to set up calculation of density derived atomic point charges.",&
            n_keywords=3, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="NUM_GAUSS",&
            description="Specifies the numbers of gaussian used to fit the QM density for each atomic site.",&
            usage="NUM_GAUSS {integer}", required=.FALSE.,&
            n_var=1, type_of_var=integer_t, default_i_val=3, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PFACTOR",&
            description="Specifies the progression factor for the gaussian exponent for each atomic site.",&
            usage="PFACTOR {real}", required=.FALSE.,&
            n_var=1, type_of_var=real_t, default_r_val=1.5_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MIN_RADIUS",&
            description="Specifies the smallest radius of the gaussian used in the fit. All other radius are"//&
            " obtained with the progression factor.",&
            usage="MIN_RADIUS {real}", required=.FALSE.,&
            n_var=1, type_of_var=real_t, default_r_val=0.5_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="RADII",&
            description="Specifies all the radius of the gaussian used in the fit for each atomic site. The use"//&
            " of this keyword disables all other keywords of this section.",&
            usage="RADII {real} {real} .. {real}", required=.FALSE.,&
            n_var=-1, type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="GCUT",&
            description="Cutoff for fit in G-space.",&
            usage="GCUT {real}", required=.FALSE.,&
            n_var=1, type_of_var=real_t, default_r_val=SQRT(6.0_dp),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_print_fitcharge_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
      
    END IF

  END SUBROUTINE create_fitcharge_section

!***************************************************************************
!!****f* input_cp2k/create_print_fitcharge_section [1.0] *
!!
!!   NAME
!!     create_print_fitcharge_section
!!
!!   FUNCTION
!!     Create the print bsse section
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
  SUBROUTINE create_print_fitcharge_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_print_fitcharge_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: print_key
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="print",&
            description="Section of possible print options in FIT_CHARGE code.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)
       
       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="Controls the printing of information regarding the run.",&
            print_level=silent_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_print_fitcharge_section

!***************************************************************************
!!****f* input_cp2k/create_bsse_section [1.0] *
!!
!!   NAME
!!     create_bsse_section
!!
!!   FUNCTION
!!     Create the BSSE section for counterpoise correction
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
  SUBROUTINE create_bsse_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_bsse_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    TYPE(section_type), POINTER :: subsection
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="BSSE",&
            description="This section is used to set up the BSSE calculation.",&
            n_keywords=3, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword, subsection)
       CALL section_create(subsection,name="FRAGMENT",&
            description="Specify the atom number belonging to this fragment.",&
            n_keywords=2, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       ! Integer
       CALL keyword_create(keyword, name="LIST",&
            description="Specifies a list of atoms.",&
            usage="LIST {integer} {integer} .. {integer}", required=.FALSE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RANGE",&
            description="Specifies the starting and the ending index of atoms range.",&
            usage="RANGE {integer} {integer}", required=.FALSE.,&
            n_var=2, type_of_var=integer_t, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_print_bsse_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
      
    END IF

  END SUBROUTINE create_bsse_section

!***************************************************************************
!!****f* input_cp2k/create_print_bsse_section [1.0] *
!!
!!   NAME
!!     create_print_bsse_section
!!
!!   FUNCTION
!!     Create the print bsse section
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
  SUBROUTINE create_print_bsse_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_print_bsse_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: print_key
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="print",&
            description="Section of possible print options in BSSE code.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)
       
       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="Controls the printing of information regarding the run.",&
            print_level=silent_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_print_bsse_section

!!****f* input_cp2k/create_metadyn_section *
!!
!!   NAME
!!     create_metadyn_section
!!
!!   FUNCTION
!!     creates the metadynamics section
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
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_motion_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_motion_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: subsection

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="motion",&
            description="This section defines a set of tool connected with the motion of the nuclei.",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection)

       CALL create_geoopt_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_md_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_metadyn_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
       
       CALL create_mc_section(subsection, error)
       call section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
       
    END IF

  END SUBROUTINE create_motion_section

!!****f* input_cp2k/create_mc_section *
!!
!!   NAME
!!     create_mc_section
!!
!!   FUNCTION
!!     creates the metadynamics section
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
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_mc_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_mc_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword    

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="mc",&
            description="This section sets parameters to set up a MonteCarlo calculation.",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)
        
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="NSTEP",&
            description="Specifies the number of MC cycles.",&
            usage="NSTEP {integer}",&
            default_i_val=100,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NMOVES",&
            description="Specifies the number of classical moves between energy evaluations. ",&
            usage="NMOVES {integer}",&
            default_i_val=4,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       

       CALL keyword_create(keyword, name="NVOLMOVES",&
            description="How many classical volume moves to bias with.",&
            usage="NVOLMOVES {integer}",&
            default_i_val=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       

       CALL keyword_create(keyword, name="NSWAPMOVES",&
            description="How many classical swap moves to try per QS swap move.",&
            usage="NSWAPMOVES {integer}",&
            default_i_val=16,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)      

       CALL keyword_create(keyword, name="LSTOP",&
            description="Makes nstep in terms of steps, instead of cycles.",&
            usage="LSTOP {logical}",&
            default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LCLASSICAL",&
            description="Bias the moves with an inexpensive classical potential.",&
            usage="LCLASSICAL {logical}",&
            default_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       

       CALL keyword_create(keyword, name="RESTART",&
            description="Read initial configuration from restart file.",&
            usage="RESTART {logical}",&
            default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       

       CALL keyword_create(keyword, name="IUPVOLUME",&
            description="Every iupvolume steps update maximum volume displacement.",&
            usage="IUPVOLUME {integer}",&
            default_i_val=10000,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       

       CALL keyword_create(keyword, name="IUPTRANS",&
            description="Every iuptrans steps update maximum translation.",&
            usage="IUPTRANS {integer}",&
            default_i_val=10000,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       

       CALL keyword_create(keyword, name="PROGRAM",&
            description="Specify the type of simulation",&
            usage="PROGRAM (TRADITIONAL|GEMC)",&
            enum_c_vals=s2a( "TRADITIONAL","GEMC"),&
            enum_i_vals=(/do_mc_traditional,do_mc_gemc/),&
            default_i_val=do_mc_traditional,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       
       
       CALL keyword_create(keyword, name="RESTART_FILE_NAME",&
            description="Name of the restart file for MC information.",&
            usage="RESTART_FILE_NAME {filename}",&
            default_c_val="XXXXXXXX",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="COORDINATE_FILE_NAME",&
            description="The file to print the current coordinates to.",&
            usage="COORDINATE_FILE_NAME {filename}",&
            default_c_val="XXXXXXXX",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ENERGY_FILE_NAME",&
            description="The file to print current energies to.",&
            usage="ENERGY_FILE_NAME {filename}",&
            default_c_val="XXXXXXXX",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DATA_FILE_NAME",&
            description="The file to print current configurational info to.",&
            usage="DATA_FILE_NAME {filename}",&
            default_c_val="XXXXXXXX",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CELL_FILE_NAME",&
            description="The file to print current cell length info to.",&
            usage="CELL_FILE_NAME {filename}",&
            default_c_val="XXXXXXXX",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_DISP_FILE_NAME",&
            description="The file to print current maximum displacement info to.",&
            usage="MAX_DISP_FILE_NAME {filename}",&
            default_c_val="XXXXXXXX",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BOX2_FILE_NAME",&
            description="For GEMC, the name of the input file for the other box.",&
            usage="BOX2_FILE_NAME {filename}",&
            default_c_val="XXXXXXXX",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PRESSURE",&
            description="The pressure for NpT simulations, in bar.",&
            usage="PRESSURE {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TEMP",&
            description="The temperature of the simulation, in Kelvin.",&
            usage="TEMP {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PMSWAP",&
            description="The probability of attempting a swap move.",&
            usage="PMSWAP {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PMVOLUME",&
            description="The probability of attempting a volume move.",&
            usage="PMVOLUME {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PMTRAION",&
            description="The probability of attempting a conformational change.",&
            usage="PMTRAION {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PMTRANS",&
            description="The probability of attempting a molecule translation.",&
            usage="PMTRANS {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMVOLUME",&
            description="Maximum volume displacement, in angstrom**3.",&
            usage="RMVOLUME {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMBOND",&
            description="Maximum bond length displacement, in angstrom.",&
            usage="RMBOND {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMANGLE",&
            description="Maximum bond angle distplacement, in degrees.",&
            usage="RMANGLE {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMROT",&
            description="Maximum rotational displacement, in degrees.",&
            usage="RMROT {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMTRANS",&
            description="Maximum translational displacement, in angstroms.",&
            usage="RMTRANS {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_mc_section

!!****f* input_cp2k/create_geoopt_section *
!!
!!   NAME
!!     create_geoopt_section
!!
!!   FUNCTION
!!     creates the geometry optimization section
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
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_geoopt_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_geoopt_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword    

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="geoopt",&
            description="This section sets the environment of the geometry optimizer.",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)
        
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="OPTIMIZER",&
            variants=(/"MINIMIZER"/),&
            description="Specify which method to use to perform a geometry optimization.",&
            usage="OPTIMIZER {BFGS|LBFGS}",&
            enum_c_vals=s2a( "BFGS","LBFGS"),&
            enum_i_vals=(/default_bfgs_method_id,default_lbfgs_method_id/),&
            default_i_val=default_bfgs_method_id,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="MAX_ITER",&
            description="Specifies the maximum number of geometry optimization steps.",&
            usage="MAX_ITER {integer}",&
            default_i_val=200,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       

       CALL keyword_create(keyword, name="MAX_H_RANK",&
            description="Maximum rank (and consequently size) of the "//&
            "approximate hessian matrix used by the LBFGS optimizer",&
            usage="MAX_H_RANK {integer}",&
            default_i_val=5,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_F_PER_ITER",&
            description="Maximum number of force evaluations per iteration"// &
            "(used for the line search)",&
            usage="MAX_F_PER_ITER {integer}",&
            default_i_val=3,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WANTED_PROJ_GRADIENT",&
            description="Requested norm thresold of the gradient multiplied "// &
            "by the approximate hessian",&
            usage="WANTED_PROJ_GRADIENT {real}",&
            default_r_val=1.0E-3_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WANTED_REL_F_ERROR",&
            description="Requested relative error on the objective function"//&
            "of the optimizer (the energy)",&
            usage="WANTED_REL_F_ERROR {real}",&
            default_r_val=1.0E-5_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_DR",&
            description="Tolerance on the maximum value of the displacement.",&
            usage="MAX_DR {real}",&
            default_r_val=0.0030_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_FORCE",&
            description="Tolerance on the maximum value of Forces.",&
            usage="MAX_FORCE {real}",&
            default_r_val=0.00045_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMS_DR",&
            description="Tolerance on RMS displacements.",&
            usage="RMS_DR {real}",&
            default_r_val=0.0015_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMS_FORCE",&
            description="Tolerance on RMS Forces.",&
            usage="RMS_FORCE {real}",&
            default_r_val=0.00030_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART",&
            description="The guess starting HESSIAN is read from the RESTART.",&
            usage="RESTART",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CONSTRAINTS",&
            variants=(/"CONSTRAINT"/),&
            description="Enables the use of constraints during the geometry optimization.",&
            usage="CONSTRAINTS",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_geoopt_section

!!****f* input_cp2k/create_metadyn_section *
!!
!!   NAME
!!     create_metadyn_section
!!
!!   FUNCTION
!!     creates the metadynamics section
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
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_metadyn_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_metadyn_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: subsection
    TYPE(keyword_type), pointer :: keyword    

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="metadyn",&
            description="This section sets parameters to set up a calculation of metadynamics.",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)
        
       NULLIFY(subsection,keyword)

       CALL keyword_create(keyword, name="NCOLVAR",&
            variants=(/"N_COLVAR"/),&
            description="Specifies the number of collective variables to use in the metadynamics."//&
            "Defaults 1 variable.",&
            usage="NCOLVAR {integer}",&
            default_i_val=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NT_HILLS",&
            description="Specify the number the interval of time step to spawn the hills.",&
            usage="NT_HILLS {integer}",&
            default_i_val=30,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TEMPERATURE",&
            variants=(/"TEMPCONTROL"/),&
            description="If a Lagrangian scheme is used the temperature control for the collective "//&
            "variables is specified. Default 100.0 K.",&
            usage="TEMPERATURE {real}",&
            default_r_val=100.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WW",&
            description="Specifies the width of the gaussian to spawn. Default 0.1 .",&
            usage="WW {real}",&
            default_r_val=0.1_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DO_HILLS",&
            description="This keyword enables the spawning of the hills. Default .FALSE.",&
            usage="DO_HILLS",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LAGRANGE",&
            description="Specifies whether an extended-lagrangian should be used. Default .FALSE.",&
            usage="LAGRANGE",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword, name="RESTART_HILLS",&
            variants=(/"RESTART"/),&
            description="Read hills from the RESTART file.",&
            usage="RESTART_HILLS",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_colvar_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_metadyn_section


!!****f* input_cp2k/create_colvar_section *
!!
!!   NAME
!!     create_colvar_section
!!
!!   FUNCTION
!!     creates the colvar section
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
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_colvar_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_colvar_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: subsection
    TYPE(keyword_type), pointer :: keyword    

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="colvar",&
            description="This section specify the nature of the collective variables.",&
            n_keywords=1, n_subsections=1, repeats=.TRUE., required=.TRUE.,&
            error=error)
        
       NULLIFY(subsection,keyword)

       CALL keyword_create(keyword, name="LAMBDA",&
            description="Specifies the lambda parameter of the collective variable in the"//&
            " extended lagrangian scheme.",&
            usage="LAMBDA {real}",&
            default_r_val=-HUGE(0.0_dp),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MASS",&
            description="Specifies the mass parameter of the collective variable in the"//&
            " extended lagrangian scheme. Default 0.0 .",&
            usage="MASS {real}",&
            default_r_val=-HUGE(0.0_dp),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SCALE",&
            description="Specifies the scale factor for the following collective variable.",&
            usage="SCALE {real}",&
            default_r_val=-HUGE(0.0_dp),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WALL_PLUS",&
            variants=(/"WALLP"/),&
            description="Activates the reflective wall on the upper limit of the collective variable.",&
            usage="WALL_PLUS {real}",&
            default_r_val=-HUGE(0.0_dp),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WALL_MINUS",&
            variants=(/"WALLM"/),&
            description="Activates the reflective wall on the lower limit of the collective variable.",&
            usage="WALL_MINUS {real}",&
            default_r_val=-HUGE(0.0_dp),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SUBSYS",&
            description="Specifies the subsystem on which to apply the collective variables. Default 1.",&
            usage="SUBSYS {integer}",&
            default_i_val=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       
    
       CALL create_colvar_var_section(subsection, section, error)
       
    END IF

  END SUBROUTINE create_colvar_section


!!****f* input_cp2k/create_colvar_var_section *
!!
!!   NAME
!!     create_colvar_var_section
!!
!!   FUNCTION
!!     creates the collective variables for the colvar section
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
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_colvar_var_section(subsection, section, error)
    TYPE(section_type), POINTER :: section, subsection
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_colvar_var_section',&
         routineP=moduleN//':'//routineN

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(subsection),cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,error,failure)

       CALL create_colvar_dist_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)      

       CALL create_colvar_angle_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)      

       CALL create_colvar_torsion_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)      

       CALL create_colvar_coord_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)      

    END IF
  END SUBROUTINE create_colvar_var_section

!!****f* input_cp2k/create_colvar_coord_section *
!!
!!   NAME
!!     create_colvar_coord_section
!!
!!   FUNCTION
!!     creates the colvar section regarded to the collective variables dist
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
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_colvar_coord_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_colvar_coord_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword    

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="coordination",&
            description="Section to define the coordination as a collective variables.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)
        
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS_FROM",&
            description="Specify indexes of atoms building the coordination variable. ",&
            usage="ATOMS_FROM {integer} {integer} ..",required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS_TO",&
            description="Specify indexes of atoms building the coordination variable. ",&
            usage="ATOMS_TO {integer} {integer} ..",required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="KINDS_FROM",&
            description="Specify kinds of atoms building the coordination variable. ",&
            usage="KINDS_FROM {integer} {integer} ..",required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="KINDS_TO",&
            description="Specify kinds of atoms building the coordination variable. ",&
            usage="KINDS_TO {integer} {integer} ..",required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="R0",&
            variants=(/"R_0"/),&
            description="Specify the R0 parameter in the coordination function.",&
            usage="R0 {real}",default_r_val=3.0_dp,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NN",&
            variants=(/"EXPON_NUMERATOR"/),&
            description="Sets the value of the numerator of the exponential factor"//&
            "in the coordination FUNCTION.",&
            usage="NN {integer}",default_i_val=6,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ND",&
            variants=(/"EXPON_DENOMINATOR"/),&
            description="Sets the value of the denominator of the exponential factor"//&
            "in the coordination FUNCTION.",&
            usage="ND {integer}",default_i_val=12,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_colvar_coord_section


!!****f* input_cp2k/create_colvar_torsion_section *
!!
!!   NAME
!!     create_colvar_torsion_section
!!
!!   FUNCTION
!!     creates the colvar section regarded to the collective variables dist
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
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_colvar_torsion_section(section,error,description)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: description
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_colvar_torsion_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword    

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       IF (PRESENT(description)) THEN
          CALL section_create(section,name="torsion",&
               description=description,&
               n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
               error=error)
       ELSE
          CALL section_create(section,name="torsion",&
               description="Section to define the torsion as a collective variables.",&
               n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
               error=error)
       END IF
        
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies the atom indexes defining the torsion.",&
            usage="ATOMS {integer} {integer} {integer} {integer}",required=.TRUE.,&
            n_var=4, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_colvar_torsion_section

!!****f* input_cp2k/create_colvar_angle_section *
!!
!!   NAME
!!     create_colvar_angle_section
!!
!!   FUNCTION
!!     creates the colvar section regarded to the collective variables dist
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
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_colvar_angle_section(section,error,description)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: description
    
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_colvar_angle_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword    

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       IF (PRESENT(description)) THEN
          CALL section_create(section,name="angle",&
               description=description,&
               n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
               error=error)
       ELSE
          CALL section_create(section,name="angle",&
               description="Section to define the angle as a collective variables.",&
               n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
               error=error)
       END IF

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies the atom indexes defining the angle.",&
            usage="ATOMS {integer} {integer} {integer}",required=.TRUE.,&
            n_var=3, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_colvar_angle_section


!!****f* input_cp2k/create_colvar_dist_section *
!!
!!   NAME
!!     create_colvar_dist_section
!!
!!   FUNCTION
!!     creates the colvar section regarded to the collective variables dist
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
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_colvar_dist_section(section,error,description)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: description

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_colvar_dist_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword    

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       IF (PRESENT(description)) THEN
          CALL section_create(section,name="distance",&
               description=description,&
               n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
               error=error)
       ELSE
          CALL section_create(section,name="distance",&
               description="Section to define the distance as a collective variables.",&
               n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
               error=error)
       END IF

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies the atom indexes defining the distance.",&
            usage="ATOMS {integer} {integer}",required=.TRUE.,&
            n_var=2, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_colvar_dist_section


!***************************************************************************
!!****f* input_cp2k/create_mm_section [1.0] *
!!
!!   NAME
!!     create_mm_section
!!
!!   FUNCTION
!!     Create the input section for FIST.. Come on.. Let's get woohooo
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
  SUBROUTINE create_mm_section(section,error)
    TYPE(section_type), pointer :: section
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_mm_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), pointer :: subsection

    failure=.false.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="mm",&
            description="This section contains all information to run a MM calculation.",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection)

       CALL create_ewald_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
       
       CALL create_forcefield_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
       
       CALL create_poisson_mm_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_print_mm_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_mm_section

!***************************************************************************
!!****f* input_cp2k/create_print_mm_section [1.0] *
!!
!!   NAME
!!     create_print_mm_section
!!
!!   FUNCTION
!!     Create the print mm section
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
  SUBROUTINE create_print_mm_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_print_mm_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: print_key
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="print",&
            description="Section of possible print options in MM code.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(print_key)

       CALL cp_print_key_section_create(print_key,"ATOMIC_COORDINATES",&
            description="Controls the printing of atomic coordinates during MM.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"BASIC_DATA_TYPES",&
            description="Controls the printing of the basic data types.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"BLACS_INFORMATION",&
            description="Controls the printing of information about the specified BLACS context.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"CELL_PARAMETERS",&
            description="Controls the printing of cell information", &
            print_level=silent_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"CPUTIME",&
            description="Controls the printing of cputime information along the MM execution.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"DERIVATIVES",&
            description="Controls the printing of derivatives.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"DIPOLE",&
            description="Controls the printing of dipole information. Computed with"//&
            " the Berry phase formalism. OutPut in atomic Units.", &
            print_level=silent_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"DISTRIBUTION",&
            description="Controls the printing of the distribution of the simple pair neighbour list.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"INTERATOMIC_DISTANCES",&
            description="Controls the printing of interatomic distances.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"KIND_INFORMATION",&
            description="Controls the printing of information on atomic kinds.", &
            print_level=low_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"MEMORY",&
            description="Controls the printing of information on memory usage", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"NEIGHBOUR_LISTS",&
            description="Controls the printing of the set of neighbor lists.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PHYSICAL_CONSTANTS",&
            description="Controls the printing of physical constants used in the code.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"POTENTIAL",&
            description="Controls the printing of the potential generated by the MM "//&
            "charge distribution on a cube file.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PRECISION",&
            description="Controls the precision of the results printed.",&
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PROGRAM_BANNER",&
            description="Controls the printing of the banner of the MM program",&
            print_level=silent_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="Controls the printing of information regarding the run.",&
            print_level=silent_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PW_GRID_INFORMATION",&
            description="Controls the printing of information regarding the pw_grid structure.",&
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"TOTAL_NUMBERS",&
            description="Controls the printing of the total number of kinds, atoms, etc.",&
            print_level=low_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_print_mm_section

!***************************************************************************
!!****f* input_cp2k/create_forcefield_section [1.0] *
!!
!!   NAME
!!     create_forcefield_section
!!
!!   FUNCTION
!!     Create the forcefield section. This section is useful to set up the
!!     proper force_field for FIST calculations
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
  SUBROUTINE create_forcefield_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_forcefield_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    TYPE(section_type), POINTER :: subsection
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="FORCEFIELD",&
            description="Section specifying information regarding how to set up properly"// &
            " a force_field for the classical calculations.",&
            n_keywords=2, n_subsections=2, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection,keyword)

       CALL keyword_create(keyword, name="PARMTYPE",&
            description="Specifies the type of parameters file to use.",&
            usage="PARMTYPE {OFF|CHM}",default_c_val="OFF",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PARMFILE",&
            description="Specifies the filename that contains the parameters of the FF.",&
            usage="PARMFILE {FILENAME}",default_c_val=" ",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RCUT_NB",&
            description="Cutoff radius for nonbonded interactions",&
            usage="RCUT_NB {real}", type_of_var=real_t, default_r_val=-1.0_dp,&
            unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)                     

       CALL keyword_create(keyword, name="VDW_SCALE14",&
            description="Scaling factor for the VDW 1-4 ",&
            usage="VDW_SCALE14 1.0", default_r_val=1.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EI_SCALE14",&
            description="Scaling factor for the electrostatics 1-4 ",&
            usage="EI_SCALE14 1.0", default_r_val=0.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       !
       ! subsections
       !
       CALL create_NONBONDED_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_NONBONDED14_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_CHARGE_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_BONDS_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_BENDS_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_TORSIONS_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_forcefield_section

!***************************************************************************
!!****f* input_cp2k/create_TORSIONS_section [1.0] *
!!
!!   NAME
!!     create_ipbv_section
!!
!!   FUNCTION
!!     This section specifies the charges of the MM atoms
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
  SUBROUTINE create_TORSIONS_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_TORSIONS_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="TORSIONS",&
            description="This section specifies the torsion potential of the MM system."//&
            "JORGENSEN type potential supported.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kinds involved in the bond.",&
            usage="ATOMS {KIND1} {KIND2} {KIND3} {KIND4}", required=.TRUE., type_of_var=char_t,&
            n_var=4, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="K",&
            description="Defines the force constant of the potential",&
            usage="K {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="hartree",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PHI0",&
            description="Defines the phase of the potential",&
            usage="PHI0 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="rad", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="M",&
            description="Defines the multiplicity of the potential",&
            usage="M {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_TORSIONS_section

!***************************************************************************
!!****f* input_cp2k/create_BENDS_section [1.0] *
!!
!!   NAME
!!     create_ipbv_section
!!
!!   FUNCTION
!!     This section specifies the charges of the MM atoms
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
  SUBROUTINE create_BENDS_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_BENDS_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="BENDS",&
            description="This section specifies the bends potential of the MM system."//&
            "HARMONIC type potential supported.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kinds involved in the bond.",&
            usage="ATOMS {KIND1} {KIND2} {KIND3}", required=.TRUE., type_of_var=char_t,&
            n_var=3, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword, name="K",&
            description="Defines the force constant of the potential",&
            usage="K {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="hartree*rad^-2",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="THETA0",&
            description="Defines the equilibrium angle.",&
            usage="THETA0 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str='rad',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_BENDS_section

!***************************************************************************
!!****f* input_cp2k/create_BONDS_section [1.0] *
!!
!!   NAME
!!     create_ipbv_section
!!
!!   FUNCTION
!!     This section specifies the charges of the MM atoms
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
  SUBROUTINE create_BONDS_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_BONDS_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="BONDS",&
            description="This section specifies the potential of the bonds of the MM system",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kinds involved in the bond.",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="K",&
            description="Defines the force constant of the potential. For QUARTIC potentials"//&
            " three numbers are expected.",&
            usage="K {real}", required=.TRUE., type_of_var=real_t,&
            n_var=-1, unit_str="hartree*bohr^-2",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="R0",&
            description="Defines the equilibrium distance.",&
            usage="R0 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="bohr", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_BONDS_section

!***************************************************************************
!!****f* input_cp2k/create_charge_section [1.0] *
!!
!!   NAME
!!     create_ipbv_section
!!
!!   FUNCTION
!!     This section specifies the charges of the MM atoms
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
  SUBROUTINE create_charge_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_charge_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="charges",&
            description="This section specifies the charges of the MM atoms",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ATOM",&
            description="Defines the atomic kind of the charge.",&
            usage="ATOM {KIND1}", required=.TRUE., type_of_var=char_t,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CHARGE",&
            description="Defines the charge of the MM atom in electron charge unit.",&
            usage="CHARGE {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_charge_section

!***************************************************************************
!!****f* input_cp2k/create_NONBONDED14_section [1.0] *
!!
!!   NAME
!!     create_NONBONDED14_section
!!
!!   FUNCTION
!!     This section specifies the input parameters for 1-4 NON-BONDED 
!!     Interactions
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
  SUBROUTINE create_NONBONDED14_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_NONBONDED14_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: subsection
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="nonbonded14",&
            description="This section specifies the input parameters for 1-4 NON-BONDED interactions.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection)
       
       CALL create_LJ_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_Williams_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_Goodwin_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_NONBONDED14_section

!***************************************************************************
!!****f* input_cp2k/create_NONBONDED_section [1.0] *
!!
!!   NAME
!!     create_NONBONDED_section
!!
!!   FUNCTION
!!     This section specifies the input parameters for 1-4 NON-BONDED 
!!     Interactions
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
  SUBROUTINE create_NONBONDED_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_NONBONDED_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: subsection
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="nonbonded",&
            description="This section specifies the input parameters for 1-4 NON-BONDED interactions.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection)
       
       CALL create_LJ_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_Williams_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_EAM_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)       

       CALL create_Goodwin_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_IPBV_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_NONBONDED_section

!***************************************************************************
!!****f* input_cp2k/create_EAM_section [1.0] *
!!
!!   NAME
!!     create_EAM_section
!!
!!   FUNCTION
!!     This section specifies the input parameters for EAM  potential type
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
  SUBROUTINE create_EAM_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_EAM_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="EAM",&
            description="This section specifies the input parameters for LENNARD-JONES potential type.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ATOM",&
            description="Defines the atomic kind involved in the nonbond potential",&
            usage="ATOM {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PARMFILE",&
            description="Specifies the filename that contains the parameters of the EAM FF.",&
            usage="PARMFILE {FILENAME}",default_c_val=" ",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_EAM_section

!***************************************************************************
!!****f* input_cp2k/create_LJ_section [1.0] *
!!
!!   NAME
!!     create_LJ_section
!!
!!   FUNCTION
!!     This section specifies the input parameters for Lennard-Jones potential type
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
  SUBROUTINE create_LJ_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_LJ_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="lennard-jones",&
            description="This section specifies the input parameters for LENNARD-JONES potential type.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ATOM",&
            description="Defines the atomic kind involved in the nonbond potential",&
            usage="ATOM {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPSILON",&
            description="Defines the EPSILON parameter of the LJ potential",&
            usage="EPSILON {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="K_e", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SIGMA",&
            description="Defines the SIGMA parameter of the LJ potential",&
            usage="SIGMA {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RCUT",&
            description="Defines the cutoff parameter of the LJ potential",&
            usage="RCUT {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_LJ_section


!***************************************************************************
!!****f* input_cp2k/create_Williams_section [1.0] *
!!
!!   NAME
!!     create_Williams_section
!!
!!   FUNCTION
!!     This section specifies the input parameters for Williams potential type
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
  SUBROUTINE create_Williams_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_Williams_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="williams",&
            description="This section specifies the input parameters for WILLIAMS potential type.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ATOM",&
            description="Defines the atomic kind involved in the nonbond potential",&
            usage="ATOM {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="A",&
            description="Defines the A parameter of the Williams potential",&
            usage="A {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="K_e", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="B",&
            description="Defines the B parameter of the Williams potential",&
            usage="B {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom^-1",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="C",&
            description="Defines the C parameter of the Williams potential",&
            usage="C {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="K_e*angstrom^6",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RCUT",&
            description="Defines the cutoff parameter of the Williams potential",&
            usage="RCUT {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_Williams_section

!***************************************************************************
!!****f* input_cp2k/create_Goodwin_section [1.0] *
!!
!!   NAME
!!     create_Goodwin_section
!!
!!   FUNCTION
!!     This section specifies the input parameters for Goodwin potential type
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
  SUBROUTINE create_Goodwin_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_Goodwin_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="goodwin",&
            description="This section specifies the input parameters for GOODWIN potential type.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ATOM",&
            description="Defines the atomic kind involved in the nonbond potential",&
            usage="ATOM {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="VR0",&
            description="Defines the VR0 parameter of the Goodwin potential",&
            usage="VR0 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="K_e",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="D",&
            description="Defines the D parameter of the Goodwin potential",&
            usage="D {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DC",&
            description="Defines the DC parameter of the Goodwin potential",&
            usage="DC {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="M",&
            description="Defines the M parameter of the Goodwin potential",&
            usage="M {real}", required=.TRUE., type_of_var=integer_t,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MC",&
            description="Defines the MC parameter of the Goodwin potential",&
            usage="MC {real}", required=.TRUE., type_of_var=integer_t,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RCUT",&
            description="Defines the cutoff parameter of the Goodwin potential",&
            usage="RCUT {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_Goodwin_section

!***************************************************************************
!!****f* input_cp2k/create_ipbv_section [1.0] *
!!
!!   NAME
!!     create_ipbv_section
!!
!!   FUNCTION
!!     This section specifies the input parameters for IPBV potential type
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
  SUBROUTINE create_ipbv_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_ipbv_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="ipbv",&
            description="This section specifies the input parameters for IPBV potential type.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ATOM",&
            description="Defines the atomic kind involved in the IPBV nonbond potential",&
            usage="ATOM {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword, name="RCUT",&
            description="Defines the cutoff parameter of the IPBV potential",&
            usage="RCUT {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_ipbv_section

!***************************************************************************
!!****f* input_cp2k/create_constraint_section [1.0] *
!!
!!   NAME
!!     create_constraint_section
!!
!!   FUNCTION
!!     Create the constraint section. This section is useful to impose
!!     constraints
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
  SUBROUTINE create_constraint_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_constraint_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: subsection
    TYPE(keyword_type), POINTER :: keyword
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="constraint",&
            description="Section specifying information regarding how to impose constraints"// &
            " on the system.",&
            n_keywords=0, n_subsections=2, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection, keyword)
       !
       ! subsections
       !

       CALL keyword_create(keyword, name="SHAKE",&
            description="Set the tolerance for the shake constraint algorithm.",&
            usage="SHAKE real",&
            default_r_val=1.0E-6_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_g3x3_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_g4x6_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_fixed_atom_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_internal_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_constraint_section

!***************************************************************************
!!****f* input_cp2k/create_internal_section [1.0] *
!!
!!   NAME
!!     create_internal_section
!!
!!   FUNCTION
!!     Create the constraint section specialized on g3x3 constraints
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
  SUBROUTINE create_internal_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_internal_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    TYPE(section_type), POINTER :: subsection
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="internals",&
            description="This section is used to constraint the internals (bonds,angles,torsions)"//&
            " variables or their combination.",&
            n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection)

       CALL create_constraint_dist_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_constraint_angle_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_constraint_tors_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_comb_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_internal_section


!!****f* input_cp2k/create_constraint_tors_section *
!!
!!   NAME
!!     create_constraint_tors_section
!!
!!   FUNCTION
!!     creates the colvar section regarded to the collective variables dist
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
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_constraint_tors_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_constraint_tors_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword    

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="torsion",&
            description="Specifies the torsion as a constrained variable.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)
        
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies the atom indexes defining the torsion.",&
            usage="ATOMS {integer} {integer} {integer} {integer}",required=.TRUE.,&
            n_var=4, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ANGLE",&
            description="Specifies the value of the angle to which constraint the torsion.",&
            usage="ANGLE {real}",required=.TRUE.,&
            n_var=1, type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_constraint_tors_section

!!****f* input_cp2k/create_constraint_angle_section *
!!
!!   NAME
!!     create_constraint_angle_section
!!
!!   FUNCTION
!!     creates the colvar section regarded to the collective variables dist
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
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_constraint_angle_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_constraint_angle_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword    

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="angle",&
            description="Specifies the angle as a constrained variable.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies the atom indexes defining the angle.",&
            usage="ATOMS {integer} {integer} {integer}",required=.TRUE.,&
            n_var=3, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ANGLE",&
            description="Specifies the value of the angle to which constraint the bending motion.",&
            usage="ANGLE {real}",required=.TRUE.,&
            n_var=1, type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_constraint_angle_section


!!****f* input_cp2k/create_constraint_dist_section *
!!
!!   NAME
!!     create_constraint_dist_section
!!
!!   FUNCTION
!!     creates the colvar section regarded to the collective variables dist
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
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_constraint_dist_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_constraint_dist_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword    

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="distance",&
            description="Specifies the bond as a constrained variable.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies the atom indexes defining the distance.",&
            usage="ATOMS {integer} {integer}",required=.TRUE.,&
            n_var=2, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DISTANCE",&
            description="Specifies the value of the constrained distance",&
            usage="DISTANCE {real}",required=.TRUE.,&
            n_var=1, type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_constraint_dist_section



!!****f* input_cp2k/create_comb_section *
!!
!!   NAME
!!     create_comb_section
!!
!!   FUNCTION
!!     creates the colvar section regarded to the collective variables dist
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
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_comb_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_comb_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: subsection
    TYPE(keyword_type), pointer :: keyword    

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
          CALL section_create(section,name="COMB",&
               description="Specify a combination of internal variables as a contrained variable.",&
               n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
               error=error)

       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="MIXING_FACTOR",&
            description="Specifies the weight of each primitive (following the order"//&
            " BONDS, ANGLES, TORSIONS) in building the combined"//&
            " internal coordinate. Specifying all bonds, angles and torsions the order within"//&
            " each group is assumed to be the order in which they appear in the input file." ,&
            usage="MIXING_FACTOR {real} {real} .. {real}", required=.TRUE.,&
            n_var=-1, type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_constraint_dist_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_constraint_angle_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_constraint_tors_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_comb_section

!***************************************************************************
!!****f* input_cp2k/create_fixed_atom_section [1.0] *
!!
!!   NAME
!!     create_fixed_atom_section
!!
!!   FUNCTION
!!     Create the constraint section specialized on g3x3 constraints
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
  SUBROUTINE create_fixed_atom_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_fixed_atom_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="fixed_atoms",&
            description="This section is used to constraint the overall atomic position (X,Y,Z).",&
            n_keywords=3, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       ! Integer
       CALL keyword_create(keyword, name="LIST",&
            description="Specifies a list of atoms to freeze.",&
            usage="LIST {integer} {integer} .. {integer}", required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RANGE",&
            description="Specifies the starting and the ending index of atoms range you wish to freeze.",&
            usage="RANGE {integer} {integer}", required=.TRUE.,&
            n_var=2, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MM_SUBSYS",&
            variants=(/"PROTEIN"/),&
            description="In a QM/MM run all  MM atoms are fixed. Only the QM fragment is allowed to move."//&
            " Default .FALSE.",&
            usage="MM_SUBSYS logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="QM_SUBSYS",&
            description="In a QM/MM run all QM atoms are fixed. Only the MM fragment is allowed to move."//&
            " Default .FALSE.",&
            usage="QM_SUBSYS logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_fixed_atom_section

!***************************************************************************
!!****f* input_cp2k/create_g3x3_section [1.0] *
!!
!!   NAME
!!     create_g3x3_section
!!
!!   FUNCTION
!!     Create the constraint section specialized on g3x3 constraints
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
  SUBROUTINE create_g3x3_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_g3x3_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="g3x3",&
            description="This section is used to set 3x3 (3 atoms and 3 distances) constraints.",&
            n_keywords=3, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       
       ! Integer
       CALL keyword_create(keyword, name="MOLECULE",&
            variants=(/"MOL"/),&
            description="Specifies the molecule number on which constraint will be applied",&
            usage="MOL {integer}", required=.TRUE., n_var=1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Atoms' index on which apply the constraint", usage="ATOMS 1 3 6",&
            n_var=-1,type_of_var=integer_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Real
       CALL keyword_create(keyword, name="DISTANCES",&
            description="The constrained distances' values.",&
            usage="DISTANCES {real} {real} {real}", type_of_var=real_t, n_var=-1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)              

    END IF

  END SUBROUTINE create_g3x3_section

!***************************************************************************
!!****f* input_cp2k/create_g4x6_section [1.0] *
!!
!!   NAME
!!     create_g4x6_section
!!
!!   FUNCTION
!!     Create the constraint section specialized on g4x6 constraints
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
  SUBROUTINE create_g4x6_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_g4x6_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="g4x6",&
            description="This section is used to set 4x6 (4 atoms and 6 distances) constraints.",&
            n_keywords=3, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       
       ! Integer
       CALL keyword_create(keyword, name="MOLECULE",&
            variants=(/"MOL"/),&
            description="Specifies the molecule number on which constraint will be applied",&
            usage="MOL {integer}", required=.TRUE., n_var=1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Atoms' index on which apply the constraint", usage="ATOMS 1 3 6 4",&
            n_var=-1,type_of_var=integer_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Real
       CALL keyword_create(keyword, name="DISTANCES",&
            description="The constrained distances' values.",&
            usage="DISTANCES {real} {real} {real} {real} {real} {real}",&
            type_of_var=real_t, n_var=-1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)              

    END IF

  END SUBROUTINE create_g4x6_section

!***************************************************************************
!!****f* input_cp2k/create_topology_section [1.0] *
!!
!!   NAME
!!     create_topology_section
!!
!!   FUNCTION
!!     Create the topology section for FIST.. and the base is running running...
!!     Contains all information regarding topology to be read in input file..
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
  SUBROUTINE create_topology_section(section,error)
    TYPE(section_type), pointer :: section
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_topology_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword

    failure=.false.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="topology",&
            description="Section specifying information regarding how to handle the topology"// &
            " for classical runs.",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       ! Logical
       CALL keyword_create(keyword, name="CHARGE_OCCUP",&
            variants=(/ "CHARGE_O" /),&
            description="Read MM charges from the OCCUP field of PDB file.",&
            usage="CHARGE_OCCUP logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CHARGE_BETA",&
            variants=(/ "CHARGE_B" /),&
            description="Read MM charges from the BETA field of PDB file.",&
            usage="CHARGE_BETA logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PARA_RES",&
            description="For a protein, each residue is now considered a molecule",&
            usage="PARA_RES logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DUMP_TOP",&
            description="This keyword enables the dumping of the topology as a PSF file.",&
            usage="DUMP_TOP logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="REORDER",&
            description="Reorder a list of atomic coordinates into order so it can be packed correctly.",&
            usage="REORDER logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MOL_CHECK",&
            description="Check molecules have the same number of atom and names.",&
            usage="MOL_CHECK logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Character
       CALL keyword_create(keyword, name="COORD_FILE",&
            variants=s2a("COORD_FI", "COORD_FILE_NAME"),&
            description="Specifies the filename that contains coordinates.",&
            usage="COORD_FILE FILENAME",required=.TRUE.,type_of_var=char_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CONN_FILE",&
            variants=(/"CONN_FIL"/),&
            description="Specifies the filename that contains the molecular connectivity.",&
            usage="CONN_FILE FILENAME",required=.TRUE.,type_of_var=char_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="COORDINATE",&
            variants=s2a( "COORDINA", "COORD"),&
            description="Set up the way in which coordinates will be read.",&
            usage="COORDINATE (PDB||CHM||OFF)", &
            enum_c_vals=s2a( "PDB","CHM","OFF"),&
            enum_i_vals=(/do_coord_pdb, do_coord_chm, do_coord_off/),&
            default_i_val=do_coord_off,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
 
       CALL keyword_create(keyword, name="CONNECTIVITY",&
            variants=(/ "CONNECTI"/),&
            description="Set up the way in which connectivities will be read. "// &
            "Default generates topology",&
            usage="CONNECTIVITY (PSF||GROMOS||GENERATE||OFF)", &
            enum_c_vals=s2a( "PSF","GROMOS","GENERATE","OFF"),&
            enum_i_vals=(/do_conn_psf, do_conn_gromos, do_conn_generate,&
                          do_conn_off/),&
            default_i_val=do_conn_generate,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXCLUDE_VDW",&
            description="Specifies which kind of Van der Waals interaction to skip.",&
            usage="EXCLUDE_VDW (1-2||1-3||1-4)", &
            enum_c_vals=s2a( "1-2", "1-3", "1-4"),&
            enum_i_vals=(/do_skip_vdw_12, do_skip_vdw_13, do_skip_vdw_14/),&
            default_i_val=do_skip_vdw_14,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXCLUDE_EI",&
            description="Specifies which kind of Electrostatic interaction to skip.",&
            usage="EXCLUDE_EI (1-2||1-3||1-4)", &
            enum_c_vals=s2a( "1-2", "1-3", "1-4"),&
            enum_i_vals=(/do_skip_ei_12, do_skip_ei_13, do_skip_ei_14/),&
            default_i_val=do_skip_ei_14,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
 
       CALL keyword_create(keyword, name="BONDPARM",&
            description="Used in conjunction with BONDPARM_FACTOR to "//&
                        "help determine wheather there is bonding "//&
                        "between two atoms based on a distance criteria. "//&
                        "Can use covalent radii information or VDW radii information",&
            usage="BONDPARM (COVALENT||VDW)", &
            enum_c_vals=s2a( "COVALENT", "VDW"),&
            enum_i_vals=(/do_bondparm_cov, do_bondparm_vdw/),&
            default_i_val=do_bondparm_cov,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NHCOPT",&
            description="Determines the number of DOF each NHC is attached to.",&
            usage="NHCOPT (GLOBAL||MOLECULE||MASSIVE||NONE)", &
            enum_c_vals=s2a( "GLOBAL", "MOLECULE", "MASSIVE", "NONE"),&
            enum_i_vals=(/do_nhcopt_global, do_nhcopt_molecule,&
                          do_nhcopt_massive, do_nhcopt_none/),&
            default_i_val=do_nhcopt_global,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       ! Integer
       CALL keyword_create(keyword, name="EXCLUDE_CHECK",&
            variants=(/"EXCLUDE_CH"/),&
            description="Used in conjunction with PARA_RES. Determines "//&
                        "the number of molecules ahead and behind that "//&
                        "needs to be checked in order to generate the "//&
                        "proper exclusion list. 3 is safe for all "//&
                        "molecules situations but 1 should work for "//&
                        "all proteins.",&
            usage="exclude_check {integer}", default_i_val=3,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Real
       CALL keyword_create(keyword, name="BONDPARM_FACTOR",&
            description="Used in conjunction with BONDPARM to help "//&
                        "determine wheather there is bonding between "//&
                        "two atoms based on a distance criteria.",&
            usage="bondparm_factor {real}", default_r_val=1.1_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_topology_section

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
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_md_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword
    TYPE(section_type), POINTER :: coeff_section, nose_section,uni_section, baro_section,&
         subsection, print_key

    failure=.false.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="md_new",&
            description="parameters needed perform an md run",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword, coeff_section, nose_section, uni_section, baro_section, subsection)
       CALL keyword_create(keyword, name="ensemble",&
            description="The ensemble that you want to sample with the md."//&
            " npt_i uses an isotropic resize of the cell, npt_f a flexible cell.",&
            usage="ensemble nve",&
            default_i_val=nve_ensemble,&
            enum_c_vals=ensemble_name_types,&
            enum_i_vals=ensemble_name_index,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="steps",&
            description="The number of steps to be performed",&
            usage="steps 100", default_i_val=3,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="timestep",&
            description="The length of an integration step.",&
            usage="timestep 1.0", default_r_val=0.5_dp,unit_str='fs',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="restart",&
            description="The way the trajectory should be done (init: from the "//&
            "initial (input file) positions, with velocities from a Boltzmann "//&
            "distribution; pos: using the latest positions stored in the restart"//&
            " file, with with velocities from a Boltzmann distribution"//&
            "; posvel: using both the latest positions and velocites stored in "//&
            "the restart file; all: restarts all.", usage="RESTART posvel", &
            default_i_val=md_init_start,&
            enum_c_vals=restart_name_types,&
            enum_i_vals=restart_name_index,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="temperature",&
            description="The temperature in K used to initialize the velocities"//&
            "with init and pos restart, and in the npt/nvt simulations", &
            usage="temperature 325.0",&
            default_r_val=300.0_dp, unit_str='K',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="temp_tol",&
            variants=(/"temp_to              ", "temperature_tolerance"/),&
            description="The maximum accepted oscillations of the temperature"//&
            "from the wanted temperature before rescaling the velocites."//&
            "If it is 0 no rescaling is performed", &
            usage="temp_tol 0.0", default_r_val=0.0_dp, unit_str='K', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="annealing",&
            description="Specifies the rescaling factor for velocities after the integration"//&
            "scheme used. Automatically enables the annealing procedure.", &
            usage="annealing 0.0", default_r_val=1.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="gamma",&
            description="Gamma parameter for the Langevin dynamics",&
            usage="gamma 0.001",&
            default_r_val=0.0_dp, unit_str='fs^-1',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="pv_ava",&
            description="Chris has to provide docs",&
            usage="pv_ava LOGICAL",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="pv_num",&
            description="Chris has to provide docs",&
            usage="pv_num LOGICAL",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_uniaxial_section(uni_section, error=error)
       CALL section_add_subsection(section, uni_section, error=error)
       CALL section_release(uni_section,error=error)       

       CALL create_barostat_section(baro_section, error=error)
       CALL section_add_subsection(section, baro_section, error=error)
       CALL section_release(baro_section,error=error)

       CALL create_coeff_section(coeff_section, error=error)
       CALL section_add_subsection(section, coeff_section, error=error)
       CALL section_release(coeff_section,error=error)

       CALL create_nose_section(nose_section, error=error)
       CALL section_add_subsection(section, nose_section, error=error)
       CALL section_release(nose_section,error=error)

       CALL create_qs_md_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)       


       CALL section_create(subsection,name="print",&
            description="Controls the printing properties during an MD run",&
            n_keywords=0, n_subsections=1, repeats=.TRUE., required=.FALSE.,&
            error=error)

       
       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"FILES",&
            description="Controls the printing of MD information on all files.", &
            print_level=silent_print_level,each=(/1/),&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SCREEN",&
            description="Controls the printing of MD information on screen.", &
            print_level=silent_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"COORDINATE",&
            description="Controls the printing of coordinates during an MD run.", &
            print_level=silent_print_level,each=(/1/),&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"VELOCITIES",&
            description="Controls the printing of velocities during an MD run.", &
            print_level=silent_print_level,each=(/1/),&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"COEFFICIENTS",&
            description="Controls the printing of coefficients during an MD run.", &
            print_level=medium_print_level,each=(/1/),&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PRESSURE",&
            description="Controls the printing of the Pressure info during an MD run.", &
            print_level=silent_print_level,each=(/1/),&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"ENERGIES",&
            description="Controls the printing of energies during an MD run.", &
            print_level=silent_print_level,each=(/1/),&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"TEMPERATURE",&
            description="Controls the printing of temperature info during an MD run.", &
            print_level=silent_print_level,each=(/1/),&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"RESTART",&
            description="Controls the dumping of the restart file during an MD run.", &
            print_level=silent_print_level,each=(/1/),&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)
       
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection, error=error)

    END IF
END SUBROUTINE create_md_section

!***************************************************************************
!!****f* input_cp2k/create_qs_md_section *
!!
!!   NAME
!!     create_qs_md_section
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: will contain the coeff section 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************

SUBROUTINE create_qs_md_section(section, error)
  IMPLICIT NONE
  TYPE(section_type), POINTER :: section
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='create_qs_md_section',&
       routineP=moduleN//':'//routineN
  TYPE(keyword_type), POINTER :: keyword

  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN  

     CALL section_create(section,name="QS_MD",&
          description="QS specific MD parameters.",&
          n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
          error=error)
       
     NULLIFY(keyword)

     CALL keyword_create(keyword, name="neighbor_lists_from_scratch",&
          variants=(/"NEIGH"/),&
          description="This keyword enables the building of the neighbouring list from scratch."//&
          " Default .TRUE.",&
          usage="NEIGH logical",&
          default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)
     
  END IF
END SUBROUTINE create_qs_md_section

!***************************************************************************
!!****f* input_cp2k/create_uniaxial_section *
!!
!!   NAME
!!     create_uniaxial_section
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: will contain the coeff section 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************

SUBROUTINE create_uniaxial_section(section, error)
  IMPLICIT NONE
  TYPE(section_type), POINTER :: section
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='create_uniaxial_section',&
       routineP=moduleN//':'//routineN
  TYPE(keyword_type), POINTER :: keyword
  TYPE(section_type), POINTER :: subsection

  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN  

       CALL section_create(section,name="uniaxial",&
            description="Parameters for uniaxial shock.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword )
       CALL keyword_create(keyword, name="PRESSURE",&
            description="Initial pressure",&
            usage="PRESSURE real",&
            type_of_var=real_t, n_var=1, unit_str='bar', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)            

        CALL keyword_create(keyword, name="ENERGY",&
            description="Initial energy",&
            usage="ENERGY real",&
            type_of_var=real_t, n_var=1, unit_str='hartree',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)         

       CALL keyword_create(keyword, name="VOLUME",&
            description="Initial compression",&
            usage="VOLUME real",&
            type_of_var=real_t, n_var=1, unit_str='angstrom^3',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)             

       CALL keyword_create(keyword, name="TEMP_B",&
            description="Temperature of the barostat.",&
            usage="TEMP_B real",&
            type_of_var=real_t, n_var=1,unit_str='K',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)             

       CALL keyword_create(keyword, name="TIMECO",&
            description="Barostat time constant",&
            usage="TIMECO real",&
            type_of_var=real_t, n_var=1,unit_str='wavenumber_t',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)             

       CALL keyword_create(keyword, name="VSHOCK",&
            description="Velocity of shock front",&
            usage="VSHOCK real",&
            type_of_var=real_t, n_var=1,unit_str='m/s',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)             

    END IF
  END SUBROUTINE create_uniaxial_section

!***************************************************************************
!!****f* input_cp2k/create_barostat_section *
!!
!!   NAME
!!     create_barostat_section
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: will contain the coeff section 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************

SUBROUTINE create_barostat_section(section, error)
  IMPLICIT NONE
  TYPE(section_type), POINTER :: section
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='create_barostat_section',&
       routineP=moduleN//':'//routineN
  TYPE(keyword_type), POINTER :: keyword

  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN  

       CALL section_create(section,name="barostat",&
            description="Parameters of barostat.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword, name="PRESSURE",&
            description="Initial pressure",&
            usage="PRESSURE real",&
            type_of_var=real_t, n_var=1,unit_str='bar',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)      

       CALL keyword_create(keyword, name="TIMECON",&
            description="Barostat time constant",&
            usage="TIMECON real",&
            type_of_var=real_t, n_var=1,unit_str='wavenumber_t',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)        

    END IF
  END SUBROUTINE create_barostat_section

!***************************************************************************
!!****f* input_cp2k/create_coeff_section *
!!
!!   NAME
!!     create_coeff_section
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: will contain the coeff section 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************

SUBROUTINE create_coeff_section(coeff_section, error)
  IMPLICIT NONE
  TYPE(section_type), POINTER :: coeff_section
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='create_coeff_section',&
       routineP=moduleN//':'//routineN
  TYPE(keyword_type), POINTER :: keyword
  TYPE(section_type), POINTER :: nose_section

  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(coeff_section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN  

       CALL section_create(coeff_section,name="kg_coeff",&
            description="coefficent parameters needed to perform a kg-pol md run",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       
       NULLIFY(nose_section, keyword)
       CALL keyword_create(keyword, name="ensemble",&
            description="md ensemble for KG-POL coefficients", usage="ENSEMBLE NVT", &
            default_i_val=nvt_ensemble,&
            enum_c_vals=(/"nve  ","nvt  "/),&
            enum_i_vals=(/nve_ensemble,nvt_ensemble/),&
            error=error)
       CALL section_add_keyword(coeff_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       call keyword_create(keyword, name="nhcopt",&
            description="global or massive thermostats for KG-POL coefficients", &
            usage="NHCOPT GLOBAL", &
            default_i_val=nhc_global,&
            enum_c_vals=s2a("global","massive"),&
            enum_desc=s2a("just one thermostat for all the coefficents",&
            "each coefficient has its own thermostat"),&
            enum_i_vals=(/nhc_global,nhc_massive/),&
            error=error)
       CALL section_add_keyword(coeff_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       call keyword_create(keyword, name="storage",&
            description="method to store the effective masses for"//&
            "KG-POL coefficients", usage="storage global", &
            default_i_val=dyn_coeff_mass_global,&
            enum_c_vals=s2a("global","kind","atom","coeff"),&
            enum_desc=s2a("just one mass for all the coefficents",&
            "each atomic kind has its own effective mass",&
            "coefficients of each atom have their own effective mass",&
            "each coefficient has its own effective mass"),&
            enum_i_vals=(/dyn_coeff_mass_global,dyn_coeff_mass_kind,&
                          dyn_coeff_mass_atom,dyn_coeff_mass_coeff/),&
            error=error)
       CALL section_add_keyword(coeff_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="temperature",&
            description="target temperature for the KG-POL coeffs",&
            default_r_val=1.0_dp, usage="temperature 1.0",&
            unit_str='K',error=error)
       CALL section_add_keyword(coeff_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="mass",&
            description="mass of the KG-POL coeffs if storage is global",&
            default_r_val=1.0_dp, usage="mass 1.0",unit_str='au_m',&
            error=error)
       CALL section_add_keyword(coeff_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="restart",&
            description="The way the coefficents should initialized.", usage="RESTART posvel", &
            default_i_val=kg_coeff_init_start,&
            enum_c_vals=s2a("init","pos","posvel"),&
            enum_desc=s2a("from scratch",&
            "position from restart, velocities from Boltzmann distribution",&
            "positions and velocities from restart"),&
            enum_i_vals=(/kg_coeff_init_start,kg_coeff_pos_restart, &
            kg_coeff_pos_vel_restart/),&
            error=error)
       CALL section_add_keyword(coeff_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_nose_section(nose_section, error=error)
       CALL section_add_subsection(coeff_section, nose_section, error=error)
       CALL section_release(nose_section,error=error)
       
    END IF
  END SUBROUTINE create_coeff_section
!***************************************************************************
!!****f* input_cp2k/create_nose_section *
!!
!!   NAME
!!     create_nose_section
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
!!     gloria
!!
!!*** **********************************************************************

SUBROUTINE create_nose_section(section, error)
  IMPLICIT NONE
  TYPE(section_type), POINTER :: section
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='create_nose_section',&
       routineP=moduleN//':'//routineN
  TYPE(keyword_type), POINTER :: keyword

  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN  
       CALL section_create(section,name="nose",&
            description="paramameters of the Nose Hoover termostat chain for KG coeffs",&
            required=.TRUE.,error=error)
       NULLIFY(keyword)
       CALL keyword_create(keyword, name="length",&
            description="length of the NoseHoover chain", usage="length integer", &
            default_i_val=2,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="Yoshida",&
            description="order of the yoshida integretor usef for the termostat",&
            usage="Yoshida integer", &
            default_i_val=3,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="timecon",&
            description="timeconstant of the termostat chain",&
            usage="timeconstant real", &
            default_r_val=1000._dp,unit_str='wavenumber_t',&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="mts", &
            variants=s2a("multiple_time_steps","mult_t_steps"),&
            description="number of multiple timesteps to be used for the NoseHoover chain",&
            usage="mts integer", &
            default_i_val=2,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

  END IF

END SUBROUTINE create_nose_section

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
         " SPME is the smooth particle mesh using beta-Euler splines."//&
         " PME is the particle mesh using fft interpolation.",&
         usage="EWALD_TYPE {EWALD|SPME|PME}",&
         default_i_val=do_ewald_ewald,&
         enum_c_vals=ewald_name_types,&
         enum_i_vals=ewald_name_index,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="alpha",&
         description="alpha parameter associated with Ewald (EWALD|PME|SPME)",&
         usage="alpha .44", default_r_val=.44_dp, unit_str='angstrom',error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="gmax",&
         description="maximum number of grid points (PME and EWALD)",&
         usage="gmax 21", default_i_val=21,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="ns_max",&
         description="number of grid points on small mesh (SPME only)",&
         usage="ns_max 12", default_i_val=12,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="o_spline",&
         description="order of the beta-Euler spline (SPME only)",&
         usage="o_spline 6", default_i_val=6,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="epsilon",&
         description="tolerance of gaussians for fft interpolation (PME only)",&
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
    TYPE(section_type), pointer :: subsection

    failure=.false.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="subsys",&
            description="a subsystem: coordinates, topology, molecules and cell",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection)
       CALL create_cell_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_coord_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_kind_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_constraint_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_topology_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_subsys_section

!***************************************************************************
!!****f* input_cp2k/create_coord_section [1.0] *
!!
!!   NAME
!!     create_coord_section
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
  SUBROUTINE create_coord_section(section,error)
    TYPE(section_type), POINTER :: section, print_key
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_coord_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword

    failure=.false.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="coord",&
            description="The coordinates for simple systems (like the QM ones)"//&
            " xyz format by default. More complex systems should be given with"//&
            " an external pdb file.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="The atoms coordinates in xyz format.",&
            usage="H 0.7 0.3 0.3", type_of_var=char_t, n_var=-2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_coord_section

!***************************************************************************
!!****f* input_cp2k/create_kind_section [1.0] *
!!
!!   NAME
!!     create_kind_section
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
  SUBROUTINE create_kind_section(section,error)
    TYPE(section_type), POINTER :: section, print_key
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_kind_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword

    failure=.false.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="kind",&
            description="The description of the kind of the atoms (mostly for QM)",&
            n_keywords=5, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="The name of the kind described in this section.",&
            usage="H", default_c_val="DEFAULT", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AUX_BASIS_SET",&
            variants=s2a("AUXILIARY_BASIS_SET", "AUX_BASIS"),&
            description="The auxliliary basis set (GTO type)",&
            usage="AUX_BASIS_SET DZVP", type_of_var=char_t, &
            n_var=1, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AUX_BASIS_NORMALIZATION",&
            variants=s2a( "AUXILIARY_BASIS_NORMALIZATION", "AUX_BASIS_NORM"),&
            description="The normalization of the auxliliary basis set",&
            usage="AUX_BASIS_NORMALIZATION NO", &
            enum_c_vals=s2a( "NONE","NO","WFN","WAVEFUNCTION","FUNCTION",&
            "DENS","DENSITY","SQUARE","UNDEFINED"),&
            enum_i_vals=(/0,0,1,1,1,2,2,2,-1/),&
            default_i_val=-1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BASIS_SET",&
            variants=s2a("ORBITAL_BASIS_SET","ORB_BASIS"),&
            description="The basis set (GTO type)",&
            usage="BASIS_SET DZVP", type_of_var=char_t, &
            n_var=1, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BASIS_NORMALIZATION",&
            variants=s2a("ORBITAL_BASIS_NORMALIZATION","ORB_BASIS_NORM"),&
            description="The normalization of the auxliliary basis set",&
            usage="AUX_BASIS_NORMALIZATION NO", &
            enum_c_vals=s2a( "NONE", "NO", "WFN", "WAVEFUNCTION",&
            "FUNCTION", "DENS", "DENSITY", "SQUARE", "UNDEFINED"),&
            enum_i_vals=(/0,0,1,1,1,2,2,2,-1/),&
            default_i_val=-1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CORE_CORRECTION",&
            description="Corrects the effective nuclear charge",&
            usage="CORE_CORRECTION 1.0", n_var=1, required=.FALSE., &
            default_r_val=0.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ELEMENT",&
            variants=(/ "ELEMENT_SYMBOL" /),&
            description="The element of the actual kind "//&
            "(if not given it is inferred from the kind name)",&
            usage="ELEMENT O", type_of_var=char_t,n_var=1,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MASS",&
            variants=s2a("ATOMIC_MASS", "ATOMIC_WEIGHT", "WEIGHT"),&
            description="The mass of the atom "//&
            "(if negative or non present it is inferred from the element symbol)",&
            usage="MASS 2.0", type_of_var=real_t,n_var=1,&
            required=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="POTENTIAL",&
            variants=(/ "POT" /),&
            description="The element of the actual kind "//&
            "(if not given it is inferred from the kind name)",&
            usage="POTENTIAL ghost", type_of_var=char_t, n_var=1,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="HARD_EXP_RADIUS",&
            description="The region where the hard density is supposed to be confined (GAPW)(in Bohr)",&
            usage="HARD_EXP_RADIUS 0.9", default_r_val=0.8_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RHO0_EXP_RADIUS",&
            description="the radius which defines the atomic region where "//&
            "the hard compensation density is confined."//&
            "should be less than HARD_EXP_RADIUS (GAPW)[Bohr]",&
            usage="RHO_EXP_RADIUS 0.9", default_r_val=0.8_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LEBEDEV_GRID",&
            description="The number of points for the angular part of "//&
            "the local grid (GAPW)",&
            usage="LEBEDEV_GRID 40", default_i_val=50,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RADIAL_GRID",&
            description="The number of points for the radial part of "//&
            "the local grid (GAPW)",&
            usage="RADIAL_GRID 70", default_i_val=50,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_kind_section

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
    TYPE(section_type), POINTER :: section, print_key
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_qmmm_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: subsection
    TYPE(keyword_type), pointer :: keyword

    failure=.false.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="qmmm",&
            description="Input for QM/MM calculations.",&
            n_keywords=6, n_subsections=3, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword, subsection )

       CALL keyword_create(keyword, name="MM_POTENTIAL_FILE_NAME",&
            description="Name of the file containing the potential expansion in gaussians.",&
            usage="MM_POTENTIAL_FILE_NAME {filename}",&
            default_c_val="MM_POTENTIAL",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
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
            "Standard value set to 1.0E-10_dp",&
            usage="eps_mm_rspace real",&
            default_r_val=1.0E-10_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)     

       CALL keyword_create(keyword, name="nocenter",&
            description="This keyword disables the automatic centering of the qm system every MD step."//&
            "It centers the system only for the first step. Default .FALSE. (i.e. the QM system is centered every step)",&
            usage="nocenter LOGICAL",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="nocenter0",&
            description="This keyword disables the automatic centering of the qm system every MD step "//&
            " even for the first step. Default .FALSE. (i.e. the QM system is centered every step)",&
            usage="nocenter LOGICAL",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="nocompatibility",&
            description="This keyword disables the compatibility of QM/MM "//&
            "potential between CPMD and CP2K implementations. The compatibility"//&
            " is achieved using an MM potential of the form: Erf[x/rc]/x + (1/rc -2/(pi^1/2*rc))*Exp[-(x/rc)^2] .",&
            usage="nocompatibility LOGICAL",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="use_geep_lib",&
            description=" This keyword enables the use of the internal GEEP library to generate the "//&
            "gaussian expansion of the MM potential. Using this keyword there's no need to provide "//&
            "the MM_POTENTIAL_FILENAME. It expects a number from 2 to 15 (the number of gaussian funtions"//&
            " to be used in the expansion.",&
            usage="use_geep_lib INTEGER",&
            default_i_val=0, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="typ_center",&
            description="This keyword specifies whether the centering is performes continuosly (NONE) "//&
            "or by unit grid steps (GRID).",&
            usage="typ_center (NONE|GRID)",&
            enum_c_vals=s2a( "NONE","GRID"),&
            enum_i_vals=(/ do_center_none, do_center_grid /),&
            default_i_val=do_center_none, error=error)
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
            description="Information about the qm kind in the qm/mm scheme",&
            n_keywords=3, n_subsections=0, repeats=.true., required=.true.,&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="The qm kind",&
            usage="O",&
            n_var=1,type_of_var=char_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MM_INDEX",&
            description="The indexes of the mm atoms that have this kind. This keyword can be repeated several times "//&
            "(useful if you have to specify many indexes).",&
            usage="MM_INDEX 1 2",&
            n_var=-1,type_of_var=integer_t,repeats=.TRUE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SUBSYS",&
            description="Specifies which subsystem the qm atom belongs to.",&
            usage="SUBSYS integer",n_var=1,type_of_var=integer_t,&
            default_i_val=1,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)    

       CALL section_create(subsection,name="MM_KIND",&
            description="Information about the mm kind in the qm/mm scheme",&
            n_keywords=2, n_subsections=0, repeats=.true., required=.false.,&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="The mm kind",&
            usage="O",&
            n_var=1,type_of_var=char_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RADIUS",&
            description="Specifies the radius of the atomic kinds [angstrom]",&
            usage="RADIUS real",n_var=1,type_of_var=real_t,&
            default_r_val=RADIUS_QMMM_DEFAULT,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CORR_RADIUS",&
            description="Specifies the correction radius of the atomic kinds [angstrom]",&
            usage="RADIUS real",n_var=1,type_of_var=real_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)    

       CALL create_qmmm_periodic_section(subsection, error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)       

       CALL create_qmmm_link_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_qmmm_interp_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_qmmm_forcefield_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)       

       CALL create_print_qmmm_section(subsection, error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_qmmm_section

!!****f* input_cp2k/create_qmmm_periodic_section *
!!
!!   NAME
!!     create_qmmm_periodic_section
!!
!!   FUNCTION
!!     Input section to set QM/MM periodic boundary conditions
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     tlaino
!!
!!*** **********************************************************************
  SUBROUTINE create_qmmm_periodic_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_qmmm_periodic_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    TYPE(section_type), POINTER :: subsection
    failure=.FALSE.
    
    NULLIFY(keyword, subsection)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="PERIODIC",&
            description="Specify parameters for QM/MM periodic boundary conditions calculations",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)    

       CALL keyword_create(keyword, name="GMAX",&
            description="Specifies the maximum value of G in the reciprocal space over which perform the Ewald sum.",&
            usage="GMAX real",n_var=1,default_r_val=1.0_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="REPLICA",&
            description="Specifies the number of replica to take into consideration for the real part of the "//&
            "calculation. Default is letting the qmmm module decide how many replica you really need.",&
            usage="REPLICA integer",n_var=1,default_i_val=-1,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       
       
       CALL create_fitcharge_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)       

    END IF
  END SUBROUTINE create_qmmm_periodic_section

!!****f* input_cp2k/create_qmmm_forcefield_section *
!!
!!   NAME
!!     create_qmmm_forcefield_section
!!
!!   FUNCTION
!!     creates the qm/mm forcefield section to override to the FF specification
!!     given in the FIST input
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     tlaino
!!
!!*** **********************************************************************
  SUBROUTINE create_qmmm_forcefield_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_qmmm_forcefield_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: subsection
    failure=.FALSE.
    
    NULLIFY(subsection )
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="FORCEFIELD",&
            description="Specify information on the QM/MM forcefield",&
            n_keywords=0, n_subsections=2, repeats=.true., required=.true.,&
            error=error)    

       CALL create_qmmm_ff_nb_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_NONBONDED14_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_qmmm_forcefield_section

!!****f* input_cp2k/create_qmmm_ff_nb_section *
!!
!!   NAME
!!     create_qmmm_ff_nb_section
!!
!!   FUNCTION
!!     creates the qm/mm forcefield section to override to the FF specification
!!     given in the FIST input - NONBONDED PART
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     tlaino
!!
!!*** **********************************************************************
  SUBROUTINE create_qmmm_ff_nb_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_qmmm_ff_nb_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER ::  subsection
    failure=.FALSE.
    
    NULLIFY(subsection )
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="NONBONDED",&
            description="Specify information on the QM/MM non-bonded forcefield",&
            n_keywords=0, n_subsections=2, repeats=.true., required=.true.,&
            error=error)    

       CALL create_LJ_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
       
       CALL create_Williams_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_Goodwin_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_qmmm_ff_nb_section

!!****f* input_cp2k/create_qmmm_link_section *
!!
!!   NAME
!!     create_qmmm_link_section
!!
!!   FUNCTION
!!     creates the qm/mm link section
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     tlaino
!!
!!*** **********************************************************************
  SUBROUTINE create_qmmm_link_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_qmmm_link_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    TYPE(section_type), POINTER :: subsection
    failure=.FALSE.
    
    NULLIFY(keyword, subsection )
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="LINK",&
            description="Specify information on the QM/MM link treatment",&
            n_keywords=7, n_subsections=2, repeats=.true., required=.true.,&
            error=error)    
       
       CALL keyword_create(keyword, name="QM_INDEX",&
            variants=(/ "QM" /),&
            description="Specifies the index of the QM atom involved in the QM/MM link",&
            usage="QM_INDEX integer",n_var=1,type_of_var=integer_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="QM_KIND",&
            description="Specifies the element of the QM capping atom involved in the QM/MM link",&
            usage="QM_KIND char",n_var=1,type_of_var=char_t,&
            default_c_val="H",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MM_INDEX",&
            variants=(/ "MM" /),&
            description="Specifies the index of the MM atom involved in the QM/MM link, Default hydrogen.",&
            usage="MM_INDEX integer",n_var=1,type_of_var=integer_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RADIUS",&
            description="Overwrite the specification of the radius only for the MM atom involved in the link."//&
            "Default is to use the same radius as for the specified type.",&
            usage="RADIUS real",n_var=1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CORR_RADIUS",&
            description="Overwrite the specification of the correction radius only for the MM atom involved in the link."//&
            "Default is to use the same correction radius as for the specified type.",&
            usage="RADIUS real",n_var=1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

        CALL keyword_create(keyword, name="LINK_TYPE",&
            variants=(/ "LINK ","TYPE ","LTYPE"/),&
            description="Specifies the method to use to treat the defined QM/MM link",&
            usage="LINK_TYPE char",&
            enum_c_vals=s2a( "IMOMM","GHO","PSEUDO"),&
            enum_i_vals=(/ do_qmmm_link_imomm, do_qmmm_link_gho, do_qmmm_link_pseudo/),&
            default_i_val=do_qmmm_link_imomm, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)      

       CALL keyword_create(keyword, name="ALPHA_IMOMM",&
            variants=s2a("ALPHA"),&
            description="Specifies the scaling factor to be used in the IMOMM QM/MM link scheme",&
            usage="ALPHA_IMOMM real",n_var=1,type_of_var=real_t,&
            default_r_val=ALPHA_IMOMM_DEFAULT,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)      

       CALL keyword_create(keyword, name="QMMM_SCALE_FACTOR",&
            variants=(/ "QMMM_CHARGE_SCALE ",&
                        "QMMM_CHARGE_FACTOR",&
                        "QMMM_SCALE_CHARGE " /),&
            description="Specifies the scaling factor for the MM charge involved in the link QM/MM."//&
            "This keyword affects only the QM/MM potential, it doesn't affect the electrostatic in "//&
            "the classical part of the code"//&
            "Default 1.0 i.e. delete the charge of the MM atom of the QM/MM link bond.",&
            usage="SCALE_FACTOR real",n_var=1,type_of_var=real_t,&
            default_r_val=CHARGE_SCALE_FACTOR,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FIST_SCALE_FACTOR",&
            variants=(/ "FIST_CHARGE_SCALE ",&
                        "FIST_CHARGE_FACTOR",&
                        "FIST_SCALE_CHARGE " /),&
            description="Specifies the scaling factor for the MM charge involved in the link QM/MM."//&
            "This keyword modifies the MM charge in FIST. The modified charge will be used then also"//&
            "for the generation of the QM/MM potential."//&
            "Default 1.0 i.e. delete the charge of the MM atom of the QM/MM link bond.",&
            usage="SCALE_FACTOR real",n_var=1,type_of_var=real_t,&
            default_r_val=CHARGE_SCALE_FACTOR,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       

       CALL keyword_create(keyword, name="SUBSYS",&
            description="Specifies which subsystem the qm/mm link atoms belong to.",&
            usage="SUBSYS integer",n_var=1,type_of_var=integer_t,&
            default_i_val=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL section_create(subsection,name="MOVE_MM_CHARGE",&
            description="Specify information to move a classical charge before the"//&
            " QM/MM energies and forces evaluation",&
            n_keywords=4, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)      

       CALL keyword_create(keyword, name="ATOM_INDEX_1",&
            variants=(/ "MM1" /),&
            description="Specifies the index of the MM atom involved in the QM/MM link to be moved",&
            usage="ATOM_INDEX_1 integer",n_var=1,type_of_var=integer_t,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOM_INDEX_2",&
            variants=(/ "MM2" /),&
            description="Specifies the index of the second atom defining the direction along which "//&
            " the atom will be moved",&
            usage="ATOM_INDEX_2 integer",n_var=1,type_of_var=integer_t,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ALPHA",&
            description="Specifies the scaling factor that defines the movement along the defined direction",&
            usage="ALPHA real",n_var=1,type_of_var=real_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)      

       CALL keyword_create(keyword, name="RADIUS",&
            description="Specifies the radius used for the QM/MM electrostatic coupling after movement",&
            usage="RADIUS real",n_var=1,type_of_var=real_t, default_r_val=0.0_dp, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CORR_RADIUS",&
            description="Specifies the correction radius used for the QM/MM electrostatic coupling after movement",&
            usage="RADIUS real",n_var=1,type_of_var=real_t, default_r_val=0.0_dp, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)


       CALL section_create(subsection,name="ADD_MM_CHARGE",&
            description="Specify information to add a classical charge before the"//&
            " QM/MM energies and forces evaluation",&
            n_keywords=5, n_subsections=0, repeats=.true., required=.true.,&
            error=error)      

       CALL keyword_create(keyword, name="ATOM_INDEX_1",&
            variants=(/ "MM1" /),&
            description="Specifies the index of the first atom defining the direction along which"//&
            " the atom will be added",&
            usage="ATOM_INDEX_1 integer",n_var=1,type_of_var=integer_t,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOM_INDEX_2",&
            variants=(/ "MM2" /),&
            description="Specifies the index of the second atom defining the direction along which "//&
            " the atom will be added",&
            usage="ATOM_INDEX_2 integer",n_var=1,type_of_var=integer_t,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ALPHA",&
            description="Specifies the scaling factor that defines the movement along the defined direction",&
            usage="ALPHA real",n_var=1,type_of_var=real_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)      

       CALL keyword_create(keyword, name="RADIUS",&
            description="Specifies the radius used for the QM/MM electrostatic coupling for the added source",&
            usage="RADIUS real",default_r_val=RADIUS_QMMM_DEFAULT,n_var=1,type_of_var=real_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CORR_RADIUS",&
            description="Specifies the correction radius used for the QM/MM electrostatic coupling for the added source",&
            usage="RADIUS real",default_r_val=RADIUS_QMMM_DEFAULT,n_var=1,type_of_var=real_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CHARGE",&
            description="Specifies the charge for the added source of QM/MM potential",&
            usage="CHARGE real",default_r_val=0.0_dp,n_var=1,type_of_var=real_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)
    END IF
  END SUBROUTINE create_qmmm_link_section

!!****f* input_cp2k/create_qmmm_interp_section *
!!
!!   NAME
!!     create_qmmm_interp_section
!!
!!   FUNCTION
!!     creates the interpolation section
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     tlaino
!!
!!*** **********************************************************************
  SUBROUTINE create_qmmm_interp_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_qmmm_interp_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    TYPE(section_type), POINTER :: print_key
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="interpolator",&
            description="kind of interpolation used between the multigrids",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.false.,&
            error=error)

       NULLIFY(keyword, print_key)

       CALL keyword_create(keyword, name="kind",&
            description="the interpolator to use",&
            usage="kind spline3",&
            default_i_val=spline3_nopbc_interp,&
            enum_c_vals=s2a("linear","spline3_nopbc","spline3"),&
            enum_i_vals=(/linear_interp,&
            spline3_nopbc_interp,spline3_pbc_interp/),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="safe_computation",&
            description="if a non unrolled calculation is to be performed in parallel",&
            usage="safe_computation OFF",&
            default_l_val=.true.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="aint_precond",&
            description="the approximate inverse to use to get the starting point"//&
            " for the linear solver of the spline3 methods",&
            usage="kind spline3",&
            default_i_val=pw_interp,&
            enum_c_vals=s2a( "copy","spl3_nopbc_aint1","spl3_nopbc_precond1",&
            "spl3_nopbc_aint2","spl3_nopbc_precond2","spl3_nopbc_precond3"),&
            enum_i_vals=(/no_precond,precond_spl3_aint, precond_spl3_1,&
            precond_spl3_aint2, precond_spl3_2, precond_spl3_3/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="precond",&
            description="The preconditioner used"//&
            " for the linear solver of the spline3 methods",&
            usage="kind spline3",&
            default_i_val=pw_interp,&
            enum_c_vals=s2a("copy","spl3_nopbc_aint1","spl3_nopbc_precond1",&
            "spl3_nopbc_aint2","spl3_nopbc_precond2","spl3_nopbc_precond3"),&
            enum_i_vals=(/no_precond,precond_spl3_aint, precond_spl3_1,&
            precond_spl3_aint2, precond_spl3_2, precond_spl3_3/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="eps_x",&
            description="accuracy on the solution for spline3 the interpolators",&
            usage="eps_x 1.e-15", default_r_val=1.e-10_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="eps_r",&
            description="accuracy on the residual for spline3 the interpolators",&
            usage="eps_r 1.e-15", default_r_val=1.e-10_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="max_iter",&
            variants=(/'maxiter'/),&
            description="the maximum number of iterations",&
            usage="max_iter 200", default_i_val=100, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       nullify(print_key)
       CALL cp_print_key_section_create(print_key,"conv_info",&
            description="if convergence information about the linear solver"//&
            " of the spline methods should be printed", &
            print_level=medium_print_level,each=(/10/),filename="__STD_OUT__",&
            add_last=(/1/),error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"spl_coeffs",&
            description="outputs a cube with the coefficents calculated for "//&
            "the spline interpolation", &
            print_level=debug_print_level,each=(/1/),&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)
    END IF
  END SUBROUTINE create_qmmm_interp_section

!***************************************************************************
!!****f* input_cp2k/create_print_qmmm_section [1.0] *
!!
!!   NAME
!!     create_print_qmmm_section
!!
!!   FUNCTION
!!     Create the print qmmm section
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
  SUBROUTINE create_print_qmmm_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_print_qmmm_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: print_key
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="print",&
            description="Section of possible print options specific of the QMMM code.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"CPUTIME",&
            description="Controls the printing of cputime information along the QMMM execution.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"fist_derivatives",&
            description="Print all derivatives after FIST calculation", &
            print_level=medium_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"MEMORY",&
            description="Controls the printing of information on memory usage", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PGF",&
            description="Controls the printing of the gaussian expansion basis set of the"//&
            " electrostatic potential", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"POTENTIAL",&
            description="Controls the printing of the QMMM  potential",&
            print_level=high_print_level,each=(/1/),filename="MM_ELPOT_QMMM.cube",&
            common_iter_levels=1,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"MM_POTENTIAL",&
            description="Controls the printing of the MM unidimensional potential on file",&
            print_level=high_print_level,each=(/1/),filename="MM_ELPOT.data",&
            common_iter_levels=1,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PROGRAM_BANNER",&
            description="Controls the printing of the banner of the MM program",&
            print_level=silent_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="Controls the printing of information regarding the run.",&
            print_level=medium_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PERIODIC_INFO",&
            description="Controls the printing of information regarding the periodic boundary condition.",&
            print_level=medium_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PW_GRID_INFORMATION",&
            description="Controls the printing of information regarding the pw_grid structure.",&
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"qmmm_charges",&
            description="Print all charges generating the QM/MM potential", &
            print_level=medium_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"qmmm_derivatives",&
            description="Print all derivatives after QM/MM calculation", &
            print_level=medium_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"qmmm_link_info",&
            description="Print all information on QM/MM links", &
            print_level=medium_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)
 
       CALL cp_print_key_section_create(print_key,"qs_derivatives",&
            description="Print QM derivatives after QS calculation", &
            print_level=medium_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_print_qmmm_section

!***************************************************************************
!!****f* input_cp2k/create_poisson_mm_section [1.0] *
!!
!!   NAME
!!     create_poisson_mm_section
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
  SUBROUTINE create_poisson_mm_section(section,error)
    TYPE(section_type), pointer :: section
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_poisson_mm_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword

    failure=.false.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="poisson_mm",&
            description="Sets up the poisson resolutor.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       nullify(keyword)
       CALL keyword_create(keyword, name="POISSON_SOLVER",&
            variants=(/"POISSON", "PSOLVER"/),& 
            description="Specify which kind of solver to use to solve the Poisson equation.",&
            usage="POISSON_SOLVER char",&
            enum_c_vals=s2a( "PERIODIC", "ANALYTIC", "MT", "HOCKNEY"),&
            enum_i_vals=(/ use_periodic, use_analytic, use_mt, use_hockney/),&
            default_i_val=use_periodic, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)    

       CALL keyword_create(keyword, name="PERIODIC",&
            description="Specify the directions on wich apply PBC ",&
            usage="PERIODIC (x|y|z|xy|xz|yz|xyz|none)",&
            enum_c_vals=s2a( "x","y","z","xy","xz","yz","xyz","none"),&
            enum_i_vals=(/ use_perd_x,  use_perd_y,   use_perd_z,&
                           use_perd_xy, use_perd_xz, use_perd_yz,&
                           use_perd_xyz, use_perd_none /),&
            default_i_val=use_perd_xyz, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)     

    END IF
  END SUBROUTINE create_poisson_mm_section
!***************************************************************************

!***************************************************************************
!!****f* input_cp2k/create_poisson_dft_section [1.0] *
!!
!!   NAME
!!     create_poisson_dft_section
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
  SUBROUTINE create_poisson_dft_section(section,error)
    TYPE(section_type), pointer :: section
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_poisson_dft_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword
    TYPE(section_type), pointer :: subsection

    failure=.false.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="poisson_dft",&
            description="Sets up the poisson resolutor.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword, subsection)
       CALL keyword_create(keyword, name="POISSON_SOLVER",&
            variants=(/"POISSON", "PSOLVER"/),& 
            description="Specify which kind of solver to use to solve the Poisson equation.",&
            usage="POISSON_SOLVER char",&
            enum_c_vals=s2a( "PERIODIC", "ANALYTIC", "MT", "HOCKNEY", "MULTIPOLE"),&
            enum_i_vals=(/ use_periodic, use_analytic, use_mt, use_hockney, use_multipole/),&
            default_i_val=use_periodic, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)    

       CALL keyword_create(keyword, name="PERIODIC",&
            description="Specify the directions on wich apply PBC ",&
            usage="PERIODIC (x|y|z|xy|xz|yz|xyz|none)",&
            enum_c_vals=s2a( "x","y","z","xy","xz","yz","xyz","none"),&
            enum_i_vals=(/ use_perd_x,  use_perd_y,   use_perd_z,&
                           use_perd_xy, use_perd_xz, use_perd_yz,&
                           use_perd_xyz, use_perd_none /),&
            default_i_val=use_perd_xyz, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)     

       CALL create_mt_poisson_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_multipole_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_poisson_dft_section
!***************************************************************************

!***************************************************************************
!!****f* input_cp2k/create_multipole_section [1.0] *
!!
!!   NAME
!!     create_multipole_section
!!
!!   FUNCTION
!!     Section to set-up parameters for decoupling using the Bloechl scheme
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
  SUBROUTINE create_multipole_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_multipole_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    TYPE(section_type), POINTER :: subsection
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="FIT_CHARGE",&
            description="This section is used to set up the decoupling of QM periodic images with "//&
            "the use of density derived atomic point charges.",&
            n_keywords=3, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="NUM_GAUSS",&
            description="Specifies the numbers of gaussian used to fit the QM density for each atomic site.",&
            usage="NUM_GAUSS {integer}", required=.FALSE.,&
            n_var=1, type_of_var=integer_t, default_i_val=3, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PFACTOR",&
            description="Specifies the progression factor for the gaussian exponent for each atomic site.",&
            usage="PFACTOR {real}", required=.FALSE.,&
            n_var=1, type_of_var=real_t, default_r_val=1.5_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MIN_RADIUS",&
            description="Specifies the smallest radius of the gaussian used in the fit. All other radius are"//&
            " obtained with the progression factor.",&
            usage="MIN_RADIUS {real}", required=.FALSE.,&
            n_var=1, type_of_var=real_t, default_r_val=0.5_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="RADII",&
            description="Specifies all the radius of the gaussian used in the fit for each atomic site. The use"//&
            " of this keyword disables all other keywords of this section.",&
            usage="RADII {real} {real} .. {real}", required=.FALSE.,&
            n_var=-1, type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="GCUT",&
            description="Cutoff for fit in G-space.",&
            usage="GCUT {real}", required=.FALSE.,&
            n_var=1, type_of_var=real_t, default_r_val=SQRT(6.0_dp),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_print_multipole_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
      
    END IF

  END SUBROUTINE create_multipole_section

!***************************************************************************
!!****f* input_cp2k/create_print_multipole_section [1.0] *
!!
!!   NAME
!!     create_print_multipole_section
!!
!!   FUNCTION
!!     Create the print multipole section
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
  SUBROUTINE create_print_multipole_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_print_multipole_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: print_key
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="print",&
            description="Section of possible print options in MULTIPOLE decoupling code.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)
       
       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="Controls the printing of information regarding the run.",&
            print_level=silent_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_print_multipole_section

!***************************************************************************
!!****f* input_cp2k/create_mt_poisson_section [1.0] *
!!
!!   NAME
!!     create_mt_poisson_section
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
  SUBROUTINE create_mt_poisson_section(section,error)
    TYPE(section_type), pointer :: section
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_mt_poisson_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword

    failure=.false.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="mt",&
            description="Sets up parameters of  Martyna-Tuckerman poisson solver",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ALPHA",&
            description="Convergence parameter ALPHA*RMIN. Default value 7.0",&
            usage="ALPHA real",&
            n_var=1,default_r_val=7.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="REL_CUTOFF",&
            description="Specify the multiplicative factor for the CUTOFF keyword in MULTI_GRID "//&
            " section. The result gives the cutoff at which the 1/r non-periodic FFT3D is evaluated."//&
            "Default is 2.0",&
            usage="REL_CUTOFF real",&
            n_var=1,default_r_val=2.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_mt_poisson_section
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
    TYPE(section_type), pointer :: subsection
    failure=.false.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"cell",&
            description="parameters needed to set up the cell",&
            n_keywords=6, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword, subsection)
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


       CALL keyword_create(keyword, name="PERIODIC",&
            description="Specify the directions on wich apply PBC ",&
            usage="PERIODIC (x|y|z|xy|xz|yz|xyz|none)",&
            enum_c_vals=s2a( "x","y","z","xy","xz","yz","xyz","none"),&
            enum_i_vals=(/ use_perd_x,  use_perd_y,   use_perd_z,&
                           use_perd_xy, use_perd_xz, use_perd_yz,&
                           use_perd_xyz, use_perd_none /),&
            default_i_val=use_perd_xyz, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)     


!Integers

       CALL keyword_create(keyword, name="SUBCELLS",&
            description="Read the subcells numbers.",&
            usage="SUBCELLS integer",&
            n_var=1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

!Subsections

       CALL create_cell_ref_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_cell_section



!!****f* input_cp2k/create_cell_ref_section *
!!
!!   NAME
!!     create_cell_ref_section
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
  SUBROUTINE create_cell_ref_section(section,error)
    TYPE(section_type), pointer :: section
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_cell_ref_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword
    failure=.false.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"cell_ref",&
            description="parameters needed to set up the reference cell",&
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

!Integers

       CALL keyword_create(keyword, name="SUBCELLS",&
            description="Read the subcells numbers.",&
            usage="SUBCELLS integer",&
            n_var=1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_cell_ref_section

!!****f* input_cp2k/create_harris_section [1.0] *
!!
!!   NAME
!!     create_harris_section
!!
!!   FUNCTION
!!     This section specifies the flags for the harris functional
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
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!*** **********************************************************************
  SUBROUTINE create_harris_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_harris_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="HARRIS",&
            description="This section specifies the flags for the calculation of the harris functional",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ACTIVATE",&
            description="Activates the harris functional.",&
            usage="ACTIVATE TRUE", required=.FALSE.,&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_harris_section

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
            usage="ngrids 1", default_i_val=1, error=error)
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
            usage="progression_factor 4", default_r_val=3._dp, error=error)
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

       CALL keyword_create(keyword, name="MULTIGRID_SET",&
            description="Activate a manual setting of the multigrids",&
            usage="MULTIGRID_SET", default_l_val=.FALSE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,name="MULTIGRID_CUTOFF",&
           variants=(/"CUTOFF_LIST"/),&
           description="List of cutoff values to set up multigrids manually",&
           usage="MULTIGRID_CUTOFF 200.0 100.0 ",&
           n_var=-1,type_of_var=real_t,repeats=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(subsection)
       CALL create_interp_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)
    END IF
  END SUBROUTINE create_mgrid_section
!***************************************************************************


!!****f* input_cp2k/create_interp_section *
!!
!!   NAME
!!     create_interp_section
!!
!!   FUNCTION
!!     creates the interpolation section
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     tlaino
!!
!!*** **********************************************************************
  SUBROUTINE create_interp_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_interp_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    TYPE(section_type), POINTER :: print_key
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="interpolator",&
            description="kind of interpolation used between the multigrids",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.false.,&
            error=error)

       NULLIFY(keyword, print_key)

       CALL keyword_create(keyword, name="kind",&
            description="the interpolator to use",&
            usage="kind spline3",&
            default_i_val=pw_interp,&
            enum_c_vals=s2a("pw","linear","spline3_nopbc","spline3"),&
            enum_i_vals=(/pw_interp,linear_interp,&
            spline3_nopbc_interp,spline3_pbc_interp/),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="aint_precond",&
            description="the approximate inverse to use to get the starting point"//&
            " for the linear solver of the spline3 methods",&
            usage="aint_precond copy",&
            default_i_val=no_precond,&
            enum_c_vals=s2a( "copy","spl3_nopbc_aint","spl3_nopbc_precond1"),&
            enum_i_vals=(/no_precond,precond_spl3_aint, precond_spl3_1/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="precond",&
            description="The preconditioner used"//&
            " for the linear solver of the spline3 methods",&
            usage="PRECOND copy",&
            default_i_val=no_precond,&
            enum_c_vals=s2a("copy","spl3_nopbc_aint","spl3_nopbc_precond1"),&
            enum_i_vals=(/no_precond,precond_spl3_aint, precond_spl3_1/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="eps_x",&
            description="accuracy on the solution for spline3 the interpolators",&
            usage="eps_x 1.e-15", default_r_val=1.e-10_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="eps_r",&
            description="accuracy on the residual for spline3 the interpolators",&
            usage="eps_r 1.e-15", default_r_val=1.e-10_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="max_iter",&
            variants=(/'maxiter'/),&
            description="the maximum number of iterations",&
            usage="max_iter 200", default_i_val=100, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       nullify(print_key)
       CALL cp_print_key_section_create(print_key,"conv_info",&
            description="if convergence information about the linear solver"//&
            " of the spline methods should be printed", &
            print_level=medium_print_level,each=(/10/),filename="__STD_OUT__",&
            add_last=(/1/),error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"spl_coeffs",&
            description="outputs a cube with the coefficents calculated for "//&
            "the spline interpolation", &
            print_level=debug_print_level,each=(/1/),&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)
    END IF
  END SUBROUTINE create_interp_section

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
            " minimum, do SD in that case ",&
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
          description="Target relative uncertainty in the location of the minimum",&
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

     CALL keyword_create(keyword, name="MULT_PRECONDITIONER",&
          description="Number of matrices constituting the preconditioner,"//&
          " to take into account large differences in the eigenvalues",&
          usage="MULT_PRECONDITIONER 2",default_i_val=1,error=error)
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

!!****f* input_cp2k/create_xas_section *
!!
!!   NAME
!!     create_xas_section
!!   FUNCTION
!!     makes the x-ray absorption section
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     03.2005 created [MI]
!!
!!   SOURCE
!!*** **********************************************************************
  SUBROUTINE create_xas_section(section,error)
    TYPE(section_type), pointer :: section
    type(cp_error_type), optional, intent(inout) :: error
  
    logical :: failure
    character(len=*), parameter :: routineN='create_xas_section',&
        routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: subsection
    TYPE(keyword_type), pointer :: keyword
    failure=.false.
  
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      CALL section_create(section,"xas",&
         description="Sets the method of choice to calculate the XAS. "//&
                     "The occupied states from  which we calculate the "//&
                     "excitation should be specified. "//&
                     "Localization of the orbitals may be useful.",&
         n_keywords=10, n_subsections=1, repeats=.FALSE., required=.false.,&
         error=error)

      nullify(keyword,subsection)
      CALL keyword_create(keyword, name="METHOD",&
           variants=(/"XAS_METHOD"/),&
           description="Method to be used to calculate XAS",&
           usage="METHOD TP_HH",&
           default_i_val=xas_none,&
           enum_c_vals=s2a( "NONE", "TP_HH", "TP_FH", "TDDFT","DSCF"),&
           enum_desc=s2a("No XAS is calculated", "Transition potential halfhole",&
                         "Transition potential fullhole", &
                         "Excitations by TDDFT","DSCF calculations for the first excited state"),&
           enum_i_vals=(/xas_none,xas_tp_hh,xas_tp_fh,xas_tddft,xas_dscf/),&
           error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="DIPOLE_FORM",&
           variants=(/"DIP_FORM"/),&
           description="Type of integral to get the oscillator strengths "//&
           "in the diipole approximation",&
           usage="DIPOLE_FORM string",&
           default_i_val=xas_dip_vel,&
           enum_c_vals=s2a( "LENGTH","VELOCITY" ),&
           enum_desc=s2a("Length form <i|e r |j>","Velocity form <i|d/dr|j>"),&
           enum_i_vals=(/xas_dip_len,xas_dip_vel/),error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)
        
      CALL keyword_create(keyword, name="MAXSTEP",&
           description="Max # of steps in the xas-scf for one excitation",&
           usage="MAXSTEP 150",&
           default_i_val=150,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)
 
      CALL keyword_create(keyword, name="CONVERGENCE",&
           variants=(/"CONV"/),&
           description="Convergence criterion for the xas-scf",&
           usage="CONVERGENCE 0.00005",&
           default_r_val=0.5D-6,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)
 
      CALL keyword_create(keyword, name="EACH_STEP",&
           description="How often  the xas spectra are computed",&
           usage="EACH_STEP 1",&
           default_i_val=1,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)
 
      CALL keyword_create(keyword, name="STATE_TYPE",&
           variants=(/"TYPE"/),&
           description="Type of the orbitas that are excited for the xas spectra calculation",&
           usage="STATE_TYPE 1S",&
           default_i_val=xas_1s_type,&
           enum_c_vals=s2a( "1S", "2S"),&
           enum_desc=s2a("1s orbitals", "2s orbitals"),&
           enum_i_vals=(/xas_1s_type,xas_2s_type/),&
           error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)
 
      CALL keyword_create(keyword, name="STATE_SEARCH",&
           description="# of states where to look for the one to be excited",&
           usage="STATE_SEARCH 1",&
           default_i_val=-1,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="FNAME",&
           description="Root of the file names where to store the spectra",&
           usage="FNAME xas_1s_o",&
           default_c_val="xas_out",error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,name="ATOMS_LU_BOUNDS",&
           variants=(/"AT_LU"/),&
           description="The lower and upper index of the atoms to be excited",&
           usage="ATOMS_LU_BOUNDS integer integer",&
           n_var=2,default_i_vals=(/0,0/), type_of_var=integer_t,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,name="ATOMS_LIST",&
           variants=(/"AT_LIST"/),&
           description="Indexes of the atoms to be excited"//&
           "This keyword can be repeated several times"//&
           "(useful if you have to specify many indexes).",&
           usage="ATOMS_LIST 1 2",&
           n_var=-1,type_of_var=integer_t,repeats=.TRUE.,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="ADDED_MOS",&
            description="Number of additional MOS added spin up only",&
            usage="ADDED_MOS {integer}", default_i_val=-1,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="NGAUSS",&
            description="Number of gto's for the expansion of the sto"//&
            "of the type given by STATE_TYPE",&
            usage="NGAUSS {integer}", default_i_val=3,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="EMISSION",&
           description="Calculate also the emission spectrum",&
           usage="EMISSION",&
           default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error) 

      CALL keyword_create(keyword, name="LOCALIZE_WFN",&
           description="Check and possibly optimize the localization of the wavefunctions",&
           usage="LOCALIZE_WFN",&
           default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error) 

      CALL create_localize_section(subsection,error)
      CALL section_add_subsection(section, subsection, error=error)
      CALL section_release(subsection,error=error)

      CALL keyword_create(keyword, name="CUBES",&
           description="Print some of the wavefunctions as cube files",&
           usage="CUBES",&
           default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error) 

      CALL keyword_create(keyword,name="CUBES_LU_BOUNDS",&
           variants=(/"CUBES_LU"/),&
           description="The lower and upper index of the states to be printed as cube files",&
           usage="CUBES_LU_BOUNDS integer integer",&
           n_var=2,default_i_vals=(/0,0/), type_of_var=integer_t,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,name="CUBES_LIST",&
           description="Indexes of the states to be printed as cube files"//&
           "This keyword can be repeated several times"//&
           "(useful if you have to specify many indexes).",&
           usage="CUBES_LIST 1 2",&
           n_var=-1,type_of_var=integer_t,repeats=.TRUE.,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)


    END IF
  END SUBROUTINE create_xas_section

!!*** **********************************************************************
!!****f* input_cp2k/create_linres_section [1.0] *
!!
!!   NAME
!!     create_linres_section
!!
!!   FUNCTION
!!     creates the input structure used to activate
!!     a linear response calculation 
!!     Available properties : none
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
!!     MI
!!
!!*** **********************************************************************
  SUBROUTINE create_linres_section(section,error)
    TYPE(section_type), pointer :: section
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_linres_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword
    TYPE(section_type), pointer :: subsection

    failure=.false.
    NULLIFY(keyword,subsection)

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="linres",&
            description="The linear response is used to calculate one of the "//&
            " following properties: chemicalshift, ...(to be implemented) ",&
            n_keywords=5, n_subsections=2, repeats=.FALSE., required=.TRUE.,&
            error=error)

      CALL keyword_create(keyword, name="EPS_SCF",&
           description="target accuracy for the scf convergence",&
           usage="EPS_SCF 1.e-6", default_r_val=1.e-5_dp,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="MAX_SCF",&
           description="Maximum number of SCF iteration to be performed for one optimization",&
           usage="MAX_SCF 200", default_i_val=50,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="OPT_METHOD",&
           variants=(/"MINIMIZER"/),&
           description="Define the algorithms used to optimize the Psi1 in the scf iteration",&
           usage="OPT_METHOD CG",&
           default_i_val=lr_diis,&
           enum_c_vals=s2a( "CG", "DIIS", "SD"),&
           enum_desc=s2a("Conjugate Gradients", "DIIS", "Steepesst Descent"),&
           enum_i_vals=(/lr_cg,lr_diis,lr_sd/),&
           error=error)
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
           default_i_val=ot_ls_none,&
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

      CALL keyword_create(keyword, name="PRECONDITIONER",&
           description="Type of preconditioner to be used with all minimization schemes."//&
                       "They differ in effectiveness, cost of construction, cost of application.",&
           usage="PRECONDITIONER FULL_SINGLE",&
           default_i_val=lr_precond_none,&
           enum_c_vals=s2a("FULL_SINGLE","FULL_KINETIC","FULL_S_INVERSE","NONE"),&
           enum_desc=s2a("Based on diagonalisation","Cholesky inversion of S and T", &
                         "Cholesky inversion of S", "skip preconditioning"),&
           enum_i_vals=(/lr_precond_single,lr_precond_kinetic,lr_precond_s_inverse,&
                         lr_precond_none/), error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="MULT_PRECONDITIONER",&
           description="Number of matrices constituting the preconditioner,"//&
           " to take into account large differences in the eigenvalues",&
           usage="MULT_PRECONDITIONER 2",default_i_val=1,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="ENERGY_GAP",&
           description="Energy gap estimate [a.u.] for preconditioning",&
           usage="ENERGY_GAP 0.1",&
           default_r_val=0.2_dp,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="LOCALIZE_PSI0",&
           description="Maximize the localization of the ground state wavefunctions",&
           usage="LOCALIZE_PSI0",&
           default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      NULLIFY(subsection)
      CALL create_localize_section(subsection,error)
      CALL section_add_subsection(section, subsection, error=error)
      CALL section_release(subsection,error=error)


       NULLIFY(subsection)
       CALL create_nmr_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_linres_section

  SUBROUTINE create_nmr_section(section,error)
    TYPE(section_type), pointer :: section
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_linres_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: subsection
    TYPE(keyword_type), pointer :: keyword
 
    failure=.false.
    NULLIFY(keyword,subsection)

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="nmr",&
            description="The chemical shift is calculated by DFPT ",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="CURRENT_DENSITY",&
            description="Calculate the current density and the "//&
            " induced magnetic field on a mesh in real space",&
            usage="CURRENT_DENSITY",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FULL",&
            description="Calculate the orbital dependent correction "//&
            " to get the full current density", usage="FULL",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


    END IF

  END SUBROUTINE create_nmr_section


END MODULE input_cp2k
