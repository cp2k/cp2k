

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!


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
  USE bibliography,                    ONLY: Frigo2005
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser,                       ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE cp_units,                        ONLY: cp_unit_set_create,&
                                             cp_unit_set_release,&
                                             cp_unit_set_type,&
                                             cp_unit_to_cp2k
  USE input_constants
  USE input_cp2k_dft,                  ONLY: create_bsse_section,&
                                             create_dft_section,&
                                             create_ep_section,&
                                             create_properties_section
  USE input_cp2k_eip,                  ONLY: create_eip_section
  USE input_cp2k_mm,                   ONLY: create_mm_section
  USE input_cp2k_motion,               ONLY: create_motion_section
  USE input_cp2k_qmmm,                 ONLY: create_qmmm_section
  USE input_cp2k_subsys,               ONLY: create_subsys_section
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_parsing,                   ONLY: section_vals_parse
  USE input_section_types,             ONLY: &
       section_add_keyword, section_add_subsection, section_create, &
       section_release, section_type, section_vals_create, section_vals_get, &
       section_vals_get_subs_vals, section_vals_type
  USE input_val_types,                 ONLY: char_t,&
                                             integer_t,&
                                             logical_t,&
                                             no_t,&
                                             real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k'
  
  PUBLIC :: create_cp2k_root_section, parsed_cp2k_input,&
       create_global_section, create_cp2k_input_reading
!!***
!****************************************************************************
CONTAINS


!!****f* input_cp2k/create_cp2k_input_reading *
!!
!!   NAME
!!     create_cp2k_input_reading
!!
!!   FUNCTION
!!     reads the cp2k input from the given filepath and returns a section_vals
!!     containing the input
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - file_path: path where the input should be read
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION create_cp2k_input_reading(file_path,para_env,error) RESULT(res)
    CHARACTER(len=*), INTENT(in)             :: file_path
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(section_vals_type), POINTER         :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'create_cp2k_input_reading', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(cp_parser_type), POINTER            :: cpparser
    TYPE(cp_unit_set_type), POINTER          :: default_units
    TYPE(section_type), POINTER              :: input_structure

    CALL timeset(routineN,"I",'',handle)
    failure=.FALSE.
    NULLIFY(res)
    NULLIFY(cpparser, input_structure, default_units)
    CALL create_cp2k_root_section(input_structure,error=error)
    CALL section_vals_create(res,input_structure, error=error)
    CALL section_release(input_structure,error=error)
    CALL parser_create(cpparser,file_name=file_path,para_env=para_env,&
         error=error)
    CALL cp_unit_set_create(default_units, "OUTPUT",error=error)
    CALL section_vals_parse(res,cpparser,root_section=.FALSE.,&
         default_units=default_units,error=error)
    CALL cp_unit_set_release(default_units,error=error)
    CALL parser_release(cpparser,error=error)
    CALL timestop(0.0_dp,handle)
  END FUNCTION create_cp2k_input_reading
!***************************************************************************

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
    TYPE(section_type), POINTER              :: root_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_cp2k_root_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section

    CALL timeset(routineN,"I",'',handle)
    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(root_section),cp_failure_level,routineP,error,failure)
    IF (.NOT.  failure) THEN
       CALL section_create(root_section,name="__ROOT__",&
            description="input file of cp2k",n_keywords=0, n_subsections=10,&
            repeats=.FALSE., required=.FALSE., error=error)
       NULLIFY(section)

       CALL create_global_section(section,error=error)
       CALL section_add_subsection(root_section,section,error=error)
       CALL section_release(section,error=error)

       CALL create_test_section(section,error=error)
       CALL section_add_subsection(root_section,section,error=error)
       CALL section_release(section,error=error)

       CALL create_debug_section(section,error=error)
       CALL section_add_subsection(root_section,section,error=error)
       CALL section_release(section,error=error)

       CALL create_motion_section(section,error=error)
       CALL section_add_subsection(root_section,section,error=error)
       CALL section_release(section,error=error)

       CALL create_multi_force_section(section,error=error)
       CALL section_add_subsection(root_section,section,error=error)
       CALL section_release(section,error=error)

       CALL create_force_eval_section(section,error=error)
       CALL section_add_subsection(root_section,section,error=error)
       CALL section_release(section,error=error)

       CALL create_farming_section(section,error=error)
       CALL section_add_subsection(root_section,section,error=error)
       CALL section_release(section,error=error)

       CALL create_ext_restart_section(section,error=error)
       CALL section_add_subsection(root_section,section,error=error)
       CALL section_release(section,error=error)

       CALL create_vib_section(section,error=error)
       CALL section_add_subsection(root_section,section,error=error)
       CALL section_release(section,error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword,name="_DEFAULT_KEYWORD_",&
            description="ignores lines outside sections",&
            usage=" ", type_of_var=no_t,error=error)
       CALL section_add_keyword(root_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE create_cp2k_root_section
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
    TYPE(section_vals_type), POINTER         :: input_file
    LOGICAL, INTENT(IN), OPTIONAL            :: check_this_section
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'parsed_cp2k_input', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_check
    TYPE(section_vals_type), POINTER         :: glob_section

    failure=.FALSE.
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
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_global_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.
    NULLIFY(print_key)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="GLOBAL",&
         description="Section with general information regarding which kind "//&
         "of simulation to perform an parameters for the whole PROGRAM",&
         n_keywords=7, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
         error=error) ! to do: required=.false.

    NULLIFY(keyword)
    CALL keyword_create(keyword, name="BLACS_GRID",&
         description="how to distribute the processors on the 2d grid needed "//&
         "by BLACS (and thus SCALAPACK)", usage="BLACS_GRID SQUARE",&
         default_i_val=blacs_square,enum_c_vals=s2a("SQUARE","ROW","COLUMN"),&
         enum_desc=s2a("Distribution by marix blocks", "Distribution by matrix rows",&
         "Distribution by matrix columns"), &
         enum_i_vals=(/blacs_square,blacs_row,blacs_col/), error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="PREFERRED_FFT_LIBRARY",&
         description="Specifies the FFT library which should be preferred. "//&
                     "If it is not available, use FFTW if this is linked in, if FFTW is not available use FFTSG. "//&
                     "Contrary to earlier CP2K versions, all libraries will result in the same grids, "//&
                     "i.e. the subset of grids which all FFT libraries can transform. "//&
                     "See EXTENDED_FFT_LENGTHS if larger FFTs or grids that more precisely match a given cutoff are needed, "//&
                     "or older results need to be reproduced. "//&
                     "FFTW is often (close to) optimal, and well tested with CP2K.",&
         usage="PREFERRED_FFT_LIBRARY FFTW",&
         citations=(/Frigo2005/),&
         default_i_val=do_fft_fftw, &
         enum_i_vals=(/do_fft_sg,do_fft_fftw,do_fft_essl,do_fft_acml,do_fft_mkl,do_fft_sci/),&
         enum_c_vals=s2a("FFTSG","FFTW","FFTESSL","FFTACML","FFTMKL","FFTSCI"),&
         enum_desc=s2a("Stefan Goedecker FFT, always available, will be used in case a FFT lib is specified and not available",&
                       "a fast portable FFT library. Recommended. This will use either FFTW2 or FFTW3, as available",&
                       "available/tuned for IBM machines",&
                       "available/tuned for AMD machines (buggy & untested)",&
                       "available/tuned for INTEL machines (untested in parallel)",&
                       "available/tuned for SGI machines (untested)"),&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="FFTW_WISDOM_FILE_NAME",&
         description="The name of the file that contains wisdom (pre-planned FFTs) for use with FFTW3. "//&
                     "Using wisdom can significantly speed up the FFTs (see the FFTW homepage for details). "//&
                     "Note that wisdom is not transferable between different computer (architectures). "//&
                     "Wisdom can be generated using the fftw-wisdom tool that is part of the fftw installation. "//&
                     "cp2k/tools/cp2k-wisdom is a script that contains some additional info, and can help "//&
                     "to generate a useful default for /etc/fftw/wisdom or particular values for a given simulation.",&
         usage="FFTW_WISDOM_FILE_NAME wisdom.dat", default_lc_val="/etc/fftw/wisdom", error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="EXTENDED_FFT_LENGTHS",&
         description="Use fft library specific values for the allows number of points in FFTs. "//&
                     "Use this to obtain, with certain fft libraries, "//&
                     "larger allowed grids, or grids that more precisely match a given cutoff. "//&
                     "IMPORTANT NOTE: in this case, the actual grids used in CP2K depends on the FFT library. "//&
                     "A change of FFT library must therefore be considered equivalent to a change of basis, "//&
                     "which implies a change of total energy. ",&
         usage="EXTENDED_FFT_LENGTHS",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="ALLTOALL_SGL",&
         description="All-to-all communication (FFT) should use single precision",&
         usage="ALLTOALL_SGL YES",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
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
         description="Which program should be run to compute forces",&
         usage="PROGRAM_NAME <STRING>",&
         enum_c_vals=program_names,&
         enum_desc=program_name_desc,&
         enum_i_vals=program_name_index,&
         default_i_val=do_cp2k, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
 
    CALL keyword_create(keyword, name="PROJECT_NAME",&
         variants=(/"PROJECT"/),&
         description="Name of the project (used to build the name of the "//&
         "trajectory, and other files generated by the program)",&
         usage="PROJECT_NAME <STRING>",&
         default_c_val="PROJECT",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="OUTPUT_FILE_NAME",&
         description="Name of the output file. "//&
         "Relevant only if automatically started (through farming for example). "//&
         "If empty uses the project name as basis for it.",&
         usage="OUTPUT_FILE_NAME {filename}",default_lc_val="", error=error)
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
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="TRACE_SYNC",&
         description="For parallel TRACEd runs: synchronize the output of all processes,"//&
                     "i.e. all output written by the parent node."//&
                     "This does not work for all possible types of CP2K runs "//&
                     "(requires that timeset/timestop is called by all processes in sync, and fails otherwise)",&
         usage="TRACE_SYNC",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword,&
                        name="SEED",&
                        description="Initial seed for the global (pseudo)random number "//&
                                    "generator to create a stream of normally Gaussian "//&
                                    "distributed random numbers.",&
                        usage="SEED 4137",&
                        default_i_val=2000,&
                        error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="SAVE_MEM",&
         description="Some sections of the input structure are deallocated when not needed,"//&
         " and reallocated only when used. This reduces the required maximum memory  ",&
         usage="SAVE_MEM",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL cp_print_key_section_create(print_key, "timings", description=&
         "Controls the printing of the timing report at the end of CP2K execution", &
         print_level=silent_print_level,filename="__STD_OUT__",&
         error=error)

    CALL keyword_create(keyword,"threshold",&
         description="Specify % of CPUTIME above which the contribution will be inserted in the"//&
         " final timing report",&
         usage="TIMINGS <real>", default_r_val=0.02_dp, error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key, "REFERENCES", description=&
         "Controls the printing of the references relevant to the calculations performed", &
         print_level=silent_print_level,filename="__STD_OUT__",&
         error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="controls the global printing",&
            print_level=silent_print_level,each=(/1/),add_last=(/1/),filename="__STD_OUT__",&
            error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key, "print", description=&
         "controls the printing of physical and mathematical constants", &
         print_level=medium_print_level,filename="__STD_OUT__",&
         error=error)

    CALL keyword_create(keyword,"BASIC_DATA_TYPES",&
         description="Controls the printing of the basic data types.",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword,"physcon",&
         description="if the printkey is active prints the physical constants",&
         default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword,"spherical_harmonics",&
         description="if the printkey is active prints the spherical harmonics",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword,"RNG_MATRICES",&
                        description="Prints the transformation matrices used by the "//&
                                    " random number generator",&
                        default_l_val=.FALSE.,&
                        lone_keyword_l_val=.TRUE.,&
                        error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword,"RNG_CHECK",&
                        description="Performs a check of the global (pseudo)random "//&
                                    "number generator (RNG) and prints the result",&
                        default_l_val=.FALSE.,&
                        lone_keyword_l_val=.TRUE.,&
                        error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword,"GLOBAL_GAUSSIAN_RNG",&
                        description="Prints the initial status of the global Gaussian "//&
                                    "(pseudo)random number stream which is mostly used for "//&
                                    "the velocity initialization",&
                        default_l_val=.FALSE.,&
                        lone_keyword_l_val=.TRUE.,&
                        error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key, "PERFORMANCE", description=&
         "Controls the of standard our performance indicators", &
         print_level=high_print_level,filename="__STD_OUT__",&
         error=error)
    CALL keyword_create(keyword,"HPM",&
         description="if the printkey is active controls the printing of the hpm info (if linked as"//&
         " external library, on IBM only)",&
         default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword,"cputime",&
         description="if the printkey is active controls the printing of the cputime at important"//&
         " points during the PROGRAM run",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword,"routines",&
         description="if the printkey is active controls the printing of the routine's name"//&
         " at important points during the PROGRAM run",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword,"memory",&
         description="if the printkey is active controls the printing of the "//&
         "memory usage at important points during the PROGRAM run",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key, "warnings", description=&
         "Controls the printing of the warnings during execution.",&
         print_level=medium_print_level,filename="__STD_OUT__",&
         error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

 END SUBROUTINE create_global_section

!!****f* input_cp2k/create_test_section [1.0] *
!!
!!   NAME
!!     
!!
!!   FUNCTION
!!     section with the tests of the libraries or external code that cp2k uses
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
  SUBROUTINE create_test_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_test_section', &
      routineP = moduleN//':'//routineN

    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    CALL section_create(section,name="test",&
         description="tests to perform",&
         n_keywords=7, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
         error=error)

    NULLIFY(keyword, print_key)
    CALL keyword_create(keyword, name="memory",&
         description="tries to allocate the given amount of memory ",&
         usage="MEMORY <REAL>",default_r_val=256.e6_dp, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="copy",&
         description="Tests the performance to copy two vectors."//&
         "The results of these tests allow to determine the size of the cache "//&
         "of the CPU. This can be used to optimize the performance of the"//&
         "FFTSG library. Tests are repeated the given number of times.",&
         usage="copy 10",default_i_val=0,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="matmul",&
         description="Tests the performance of different kinds of matrix matrix "//&
         "multiply kernels for the BLAS and F95 INTRINSIC matmul. Tests are repeated the given number of times.",&
         usage="matmul 10",default_i_val=0,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="fft",&
         description="Tests the performance of all available FFT libraries for "//&
         "3D FFTs Tests are repeated the given number of times.",&
         usage="fft 10",default_i_val=0,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="pw_fft",&
         description="Tests the performance of PW calls to FFT routines. Tests are repeated the given number of times.",&
         usage="pw_fft 10",default_i_val=0,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="clebsch_gordon",variants=(/"clebsch"/),&
         description="Tests the Clebsch-Gordon Coefficients. "//&
         "Tests are repeated the given number of times.",&
         usage="clebsch_gordon 10",default_i_val=0,error=error)
    
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="rs_grids",&
         description="Tests the realspace grids. Tests are repeated the given number of times.",&
         usage="rs_grids 10",default_i_val=0,error=error)
    
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="mpi",&
         description="Tests mpi, quickly adapted benchmark code,"//&
         "will ONLY work on an even number of CPUs. comm is the relevant, "//&
         "initialized communicator. This test will produce messages "//&
         "of the size 8*10**requested_size, where requested_size is the value "//&
         "given to this keyword",&
         usage="mpi 6",default_i_val=0,error=error)
    
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="random_number_generator",variants=(/"rng"/),&
         description=" Tests the parallel random number generator (RNG)",&
         usage="rng 1000000",default_i_val=0,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="eigensolver",&
         description=" Tests the eigensolver,ABS(n) is the matrix dimension: "//&
         "n > 0: the matrix will be filled WITH random numbers; "//&
         "n < 0: the matrix is READ from the formatted file named MATRIX; "//&
         "m > 0: eigenvalues will be calculated",&
         usage="eigensolver n m",default_i_vals=(/0,0/),error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL cp_print_key_section_create(print_key,"GRID_INFORMATION",&
         description="Controls the printing of information regarding the PW and RS grid structures"//&
         " (ONLY for TEST run).",&
         print_level=medium_print_level,each=(/1/),filename="__STD_OUT__",&
         error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)    

  END SUBROUTINE create_test_section
!***************************************************************************


!!****f* input_cp2k/create_debug_section [1.0] *
!!
!!   NAME
!!     
!!
!!   FUNCTION
!!     section to setup debugging parameter
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
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_debug_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_debug_section', &
      routineP = moduleN//':'//routineN

    TYPE(keyword_type), POINTER              :: keyword

    CALL section_create(section,name="debug",&
         description="Section to setup parameters for debugging runs",&
         n_keywords=7, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
         error=error)

    NULLIFY(keyword)
    CALL keyword_create(keyword, name="numerical_pressure",&
         description="Computes the numerical tensor",&
         usage="NUMERICAL_PRESSURE <LOGICAL>",default_l_val=.FALSE.,&
         lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="DEBUG_PRESSURE",&
         description="To skip the computation of the PRESSURE in a DEBUG run",&
         usage="DEBUG_PRESSURE <LOGICAL>",default_l_val=.TRUE.,&
         lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)    

    CALL keyword_create(keyword, name="dx",&
         description=" Setup the DX parameter to evaluate numerical derivatives",&
         usage="DX <REAL>",default_r_val=0.001_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

  END SUBROUTINE create_debug_section
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
  SUBROUTINE create_multi_force_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_multi_force_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="MULTIPLE_FORCE_EVALS",&
         description="Describes how to handle multiple force_evals.",&
         n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
         error=error)

    NULLIFY(keyword)
    CALL keyword_create(keyword, name="FORCE_EVAL_ORDER",&
         description='Specify the orders of the different force_eval. When using a MIXED force_eval'//&
         " this does not need to be specified in this list, because it that takes into account only the real"//&
         " energy contributions",&
         usage="FORCE_EVAL_ORDER <INTEGER> .. <INTEGER>", type_of_var=integer_t, n_var=-1,&
         default_i_vals=(/1/), error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MULTIPLE_SUBSYS",&
         description="Specify if force_eval have different subsys. In case hey share the same subsys,"//&
         " it needs to be specified only in the MIXED force_eval (if using MIXED) or"//&
         " in the force_eval corresponding to first force_eval of the previous order (when not using MIXED).",&
         default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

   END SUBROUTINE create_multi_force_section
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
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_force_eval_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="force_eval",&
         description="parameters needed to calculate energy and forces and"//&
         " describe the system you want to analyze.",&
         n_keywords=1, n_subsections=9, repeats=.TRUE., required=.TRUE.,&
         error=error)

    NULLIFY(subsection)
    NULLIFY(keyword)
    CALL keyword_create(keyword, name="METHOD",&
         description="Which method should be used to compute forces",&
         usage="METHOD <STRING>",&
         enum_c_vals=method_names,&
         enum_desc=method_name_desc,&
         enum_i_vals=method_name_index,&
         default_i_val=do_qs, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL create_rescale_force_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_ep_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_dft_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_eip_section(subsection, error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection, error=error)

    CALL create_mm_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_mix_section(subsection,error)
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

    CALL create_f_env_print_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

  END SUBROUTINE create_force_eval_section
!***************************************************************************

!!****f* input_cp2k/create_ext_restart_section [1.0] *
!!
!!   NAME
!!     create_ext_restart_section
!!
!!   FUNCTION
!!     Creates the exteranal restart section
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
  SUBROUTINE create_ext_restart_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ext_restart_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="EXT_RESTART",&
            description="Section for external restart, specifies an external "//&
            "input file where to take positions,...",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="RESTART_FILE_NAME",variants=(/"EXTERNAL_FILE"/),&
            description="Specifies the name of restart file (or any other input file)"//&
                        " to be read. Only fields relevant to a restart will be used"//&
                        " (unless switched off with the keywords in this section)",&
            default_lc_val=" ", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART_DEFAULT",&
            description="This keyword controls the default value for all possible "//&
            " restartable keywords, unless explicitly defined. For example setting"//&
            " this keyword to .FALSE. does not restart any quantity. If, at the "//&
            " same time, one keyword is set to .TRUE. only that quantity will be"//&
            " restarted.", default_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_COUNTERS",&
            description="Restarts the counters in MD schemes",&
            type_of_var=logical_t, lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_POS",&
            description="Takes the positions from the external file",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_VEL",&
            description="Takes the velocities from the external file",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_RANDOMG",&
            description="Restarts the random number generator from the external file",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_KG_COEFF_POS",&
            description="Restarts the positions of coefficents of kg from the external file",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_KG_COEFF_VEL",&
            description="Restarts the velocity of coefficents of kg from the external file",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_SHELL_POS",&
            description="Takes the positions of the shells from the external file (only if shel-model)",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_CORE_POS",&
            description="Takes the positions of the cores from the external file (only if shel-model)",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART_SHELL_VELOCITY",&
            description="Takes the velocities of the shells from the external file (only if shel-model)",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_CORE_VELOCITY",&
            description="Takes the velocities of the shells from the external file (only if shel-model)",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_KG_THERMOSTAT",&
            description="Restarts the nose thermostats of the kg coefficents "//&
            "from the EXTERNAL file",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_BAROSTAT",&
            description="Restarts the barostat from the external file",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_BAROSTAT_THERMOSTAT",&
            description="Restarts the barostat thermostat from the external file",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_SHELL_THERMOSTAT",&
            description="Restarts the shell thermostat from the external file",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_THERMOSTAT",&
            description="Restarts the nose thermostats of the particles "//&
            "from the EXTERNAL file",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_CELL",&
            description="Restarts the cell (and cell_ref) "//&
            "from the EXTERNAL file",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_METADYNAMICS",&
            description="Restarts hills from a previous metadynamics run "//&
            "from the EXTERNAL file",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_BAND",&
            description="Restarts positions and velocities of the Band.",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_QMMM",&
            description="Restarts the following specific QMMM info: translation vectors.",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_CONSTRAINT",&
            description="Restarts constraint section. It's necessary when doing restraint"//&
            " calculation to have a perfect energy conservation. For constraints only it's"//&
            " use is optional.",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="CUSTOM_PATH",&
            description="Takes the given path "//&
            "from the EXTERNAL file",&
            type_of_var=char_t,repeats=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_ext_restart_section
!***************************************************************************

!!****f* input_cp2k/create_vib_section [1.0] *
!!
!!   NAME
!!     create_vib_section
!!
!!   FUNCTION
!!     Creates the exteranal restart section
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
  SUBROUTINE create_vib_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_vib_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="VIBRATIONAL_ANALYSIS",&
            description="Section to setup parameters to perform a Normal Modes analysis.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)
       NULLIFY(keyword, print_key)
       
       CALL keyword_create(keyword, name="DX",&
            description="Specify the increment to be used to construct the HESSIAN with "//&
            "finite difference method",&
            default_r_val=1.0E-2_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NPROC_REP",&
            description="Specify the number of processors to be used per replica "//&
            "environment (for parallel runs)",&
            default_i_val=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FULLY_PERIODIC",&
            description="Avoids to clean rotations from the Hessian matrix.",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL cp_print_key_section_create(print_key,"program_run_info",&
            description="Controls the printing basic info about the vibrational method", &
            print_level=medium_print_level,each=(/1/),add_last=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)
       
    END IF
  END SUBROUTINE create_vib_section
!***************************************************************************

!!****f* input_cp2k/create_rescale_force_section [1.0] *
!!
!!   NAME
!!     create_rescale_force_section
!!
!!   FUNCTION
!!     Creates the section controlling the rescaling of forces
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
  SUBROUTINE create_rescale_force_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_rescale_force_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="RESCALE_FORCES",&
            description="Section controlling the rescaling of forces. Useful when"//&
            " starting from quite bad geometries with unphysically large forces.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="MAX_FORCE",&
            description="Specify the Maximum Values of the force. If the force"//&
            " of one atom exceed this value it's rescaled to the MAX_FORCE"//&
            " value.",&
            default_r_val=cp_unit_to_cp2k(value=50.0_dp,&
                                          unit_str="kcalmol*angstrom^-1",&
                                          error=error),&
             unit_str="hartree*bohr^-1",&
             error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
    END IF
  END SUBROUTINE create_rescale_force_section
!***************************************************************************

!!****f* input_cp2k/create_farming_section *
!!
!!   NAME
!!     create_farming_section
!!
!!   FUNCTION
!!     creates the farming section
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
SUBROUTINE create_farming_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_farming_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: sub_section

  failure=.FALSE.
  
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  CALL section_create(section,name="farming",&
         description="describes a farming job, in which multiple inputs are executed",&
         repeats=.FALSE., required=.TRUE., error=error)
  NULLIFY(keyword)

  CALL keyword_create(keyword, name="MASTER_SLAVE",&
         description="If a master-slave setup should be employed, in which one process is used to distribute the tasks. "//&
                     "This is most useful to load-balance if not all jobs have the same length, "//&
                     "and a lot of CPUs/groups are availabe.",&
         usage="MASTER_SLAVE",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
  CALL section_add_keyword(section,keyword,error=error)
  CALL keyword_release(keyword,error=error)

  CALL keyword_create(keyword, name="NGROUPS",variants=(/"NGROUP"/),&
         description="Gives the preferred number of working groups.",&
         usage="ngroups 4", type_of_var=integer_t, error=error)
  CALL section_add_keyword(section,keyword,error=error)
  CALL keyword_release(keyword,error=error)

  CALL keyword_create(keyword, name="GROUP_SIZE",&
         description="Gives the preferred size of a working group, "//&
                     "groups will always be equal or larger than this size.",&
         usage="group_size 2", default_i_val=8, error=error)
  CALL section_add_keyword(section,keyword,error=error)
  CALL keyword_release(keyword,error=error)

  CALL keyword_create(keyword, name="GROUP_PARTITION",&
         description="gives the exact number of processors for each group.",& 
         usage="group_partition  2 2 4 2 4 ", type_of_var=integer_t, n_var=-1, error=error)
  CALL section_add_keyword(section,keyword,error=error)
  CALL keyword_release(keyword,error=error)

  CALL keyword_create(keyword, name="RESTART",&
         description="restart a farming job (and should pick up where the previous left off)",&
         usage="RESTART",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,&
         error=error)
  CALL section_add_keyword(section,keyword,error=error)
  CALL keyword_release(keyword,error=error)

  CALL keyword_create(keyword, name="MAX_JOBS_PER_GROUP",&
         variants=(/"MAX_JOBS"/),&
         description="maximum number of jobs executed per group",&
         usage="max_step 4", default_i_val=HUGE(0)/32768, error=error)
  CALL section_add_keyword(section,keyword,error=error)
  CALL keyword_release(keyword,error=error)

  CALL keyword_create(keyword, name="CYCLE",&
         description="If farming should process all jobs in a cyclic way, stopping only if MAX_JOBS_PER_GROUP is exceeded." ,&
         usage="CYCLE",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,&
         error=error)
  CALL section_add_keyword(section,keyword,error=error)
  CALL keyword_release(keyword,error=error)

  CALL keyword_create(keyword, name="WAIT_TIME",&
         description="Time to wait [s] for a new task if no task is currently available, make this zero if no clock is available",&
         usage="WAIT_TIME 0.1",default_r_val=0.5_dp, error=error)
  CALL section_add_keyword(section,keyword,error=error)
  CALL keyword_release(keyword,error=error)

  NULLIFY(sub_section)
  CALL section_create(sub_section,name="JOB",&
         description="description of the jobs to be executed",&
         repeats=.TRUE., required=.TRUE., error=error)

  CALL keyword_create(keyword,name="DIRECTORY",&
       description="the directory in which the job should be executed",&
       usage="DIRECTORY /my/path",&
       required=.TRUE.,default_lc_val=".",error=error)
  CALL section_add_keyword(sub_section,keyword,error=error)
  CALL keyword_release(keyword,error=error)

  CALL keyword_create(keyword,name="INPUT_FILE_NAME",&
       description="the filename of the input file",&
       usage="input_file my_input.inp",&
       required=.TRUE.,default_lc_val="input.inp",error=error)
  CALL section_add_keyword(sub_section,keyword,error=error)
  CALL keyword_release(keyword,error=error)
  
  CALL keyword_create(keyword,name="JOB_ID",&
       description="An ID used to indentify a job in DEPENDENCIES. "//&
                   "JOB_IDs do not need to be unique, dependencies will be on all jobs with a given ID. "//&
                   "If no JOB_ID is given, the index of the &JOB section in the input file will be used. ",&
       usage="JOB_ID 13", type_of_var=integer_t, error=error)
  CALL section_add_keyword(sub_section,keyword,error=error)
  CALL keyword_release(keyword,error=error)

  CALL keyword_create(keyword,name="DEPENDENCIES",&
       description="specifies a list of JOB_IDs on which the current job depends. "//&
                   "The current job will not be executed before all the dependencies have finished. "//&
                   "The keyword requires a MASTER_SLAVE farming run. "//&
                   "Beyond the default case, some special cases might arise: "//&
                   "1) circular dependencies will lead to a deadlock. "//&
                   "2) This keyword is not compatible with CYCLE. "//&
                   "3) MAX_JOBS_PER_GROUP is ignored (though only a total of MAX_JOBS_PER_GROUP*NGROUPS jobs will be executed) "//&
                   "4) dependencies on jobs that will not be executed (due to RESTART or MAX_JOBS_PER_GROUP) are ignored. "//&
                   "Additionally, note that, on some file systems, "//&
                   " output (restart) files might not be immediately available on all compute nodes,"//&
                   "potentially resulting in unexpected failures.", &
       usage="DEPENDENCIES 13 1 7",type_of_var=integer_t, n_var=-1, error=error)
  CALL section_add_keyword(sub_section,keyword,error=error)
  CALL keyword_release(keyword,error=error)

  CALL section_add_subsection(section, sub_section, error=error)
  CALL section_release(sub_section,error=error)

END SUBROUTINE create_farming_section
!***************************************************************************

!!****f* input_cp2k/create_f_env_print_section *
!!
!!   NAME
!!     create_f_env_print_section
!!
!!   FUNCTION
!!     -
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
SUBROUTINE create_f_env_print_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_f_env_print_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: print_key

  failure=.FALSE.

  NULLIFY(print_key)
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL section_create(section,name="print",&
          description="propreties that you want to output and that are common to all methods",&
          n_keywords=0, n_subsections=5, repeats=.FALSE., required=.TRUE.,&
          error=error)

     CALL cp_print_key_section_create(print_key,"FORCES",&
          description="Prints the atomic forces at the end of each force_eval",&
          print_level=high_print_level, filename="__STD_OUT__",error=error)
     CALL section_add_subsection(section,print_key,error=error)
     CALL section_release(print_key,error=error)

     CALL cp_print_key_section_create(print_key,"GRID_INFORMATION",&
          description="Controls the printing of information regarding the PW and RS grid structures.",&
          print_level=medium_print_level,each=(/1/),filename="__STD_OUT__",&
          error=error)
     CALL section_add_subsection(section,print_key,error=error)
     CALL section_release(print_key,error=error)

     CALL cp_print_key_section_create(print_key,"TOTAL_NUMBERS",&
          description="controls the printing of the total number of atoms, kinds,...",&
          print_level=medium_print_level, filename="__STD_OUT__",error=error)
     CALL section_add_subsection(section,print_key,error=error)
     CALL section_release(print_key,error=error)

     CALL cp_print_key_section_create(print_key,"DISTRIBUTION",&
          description="controls the printing of the distribution of molecules, atoms, ...",&
          print_level=medium_print_level, filename="__STD_OUT__",error=error)
     CALL section_add_subsection(section,print_key,error=error)
     CALL section_release(print_key,error=error)

     CALL cp_print_key_section_create(print_key,"DISTRIBUTION2D",&
          description="controls the printing of the distribution of matrix blocks,...",&
          print_level=high_print_level, filename="__STD_OUT__",error=error)
     CALL section_add_subsection(section,print_key,error=error)
     CALL section_release(print_key,error=error)

     CALL cp_print_key_section_create(print_key,"DISTRIBUTION1D",&
          description="Each node prints out its distribution info...",&
          print_level=high_print_level, filename="__STD_OUT__",error=error)
     CALL section_add_subsection(section,print_key,error=error)
     CALL section_release(print_key,error=error)

  END IF
END SUBROUTINE create_f_env_print_section
!***************************************************************************

!!****f* input_cp2k/create_cpot_section *
!!
!!   NAME
!!     create_cpot_section
!!
!!   FUNCTION
!!     section for the chemical potential difference run
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: will contain the cpot section
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE create_cpot_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_cpot_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

  failure=.FALSE.

  NULLIFY(keyword)
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL section_create(section,name="cpot",&
          description="controls a cpot (chemical potential) run",&
          n_keywords=0, n_subsections=5, repeats=.FALSE., required=.TRUE.,&
          error=error)

     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)
     CALL keyword_create(keyword,name="state_1_2_input",&
          description="the input file for the second box of the first state",&
          usage="state_1_2_input my_input.inp",&
          required=.TRUE.,default_lc_val="s-1-2.inp",error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)
     CALL keyword_create(keyword,name="state_2_1_input",&
          description="input file for the box 1 of state 2",&
          usage="state_2_1_input my_input.inp",&
          required=.TRUE.,default_lc_val="s-1-1.inp",error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)
     CALL keyword_create(keyword,name="state_2_2_input",&
          description="the input file for the second box of the second state",&
          usage="state_2_2_input my_input.inp",&
          required=.TRUE.,default_lc_val="s-2-2.inp",error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)
     CALL keyword_create(keyword,name="statistics_after",&
          description="number of steps after which the collection of statistics begins",&
          usage="statistics_after 1000",&
          required=.FALSE.,default_i_val=0,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)
     CALL keyword_create(keyword,name="statistics_each",&
          description="how often statistic is collected",&
          usage="statistics_each 5",&
          required=.FALSE.,default_i_val=0,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)
     CALL keyword_create(keyword, name="TEMP",&
          description="The temperature of the simulation, in Kelvin.",&
          usage="TEMP {real}",&
          type_of_var=real_t, required=.TRUE.,unit_str="K", error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)
  END IF
END SUBROUTINE create_cpot_section
!***************************************************************************
!***************************************************************************

!!***create_mix_section [1.0] *
!!
!!   NAME
!!     create_mixed_section
!!
!!   FUNCTION
!!     Create the input section for MIXED.
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
!!     fschiff
!!
!!*** **********************************************************************
  SUBROUTINE create_mix_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mix_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, sub2section, &
                                                sub3section, subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="MIXED",&
            description="This section contains all information to run with a hamiltonian "//&
                        "defined by a mixing of force_evals",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword, subsection, print_key)

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
            enum_i_vals=(/mix_linear_combination,mix_minimum,mix_coupled,mix_restrained,mix_generic/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="GROUP_PARTITION",&
            description="gives the exact number of processors for each group."//&
            " If not specified processors allocated will be equally distributed for"//&
            " the specified subforce_eval, trying to build a number of groups equal to the"//&
            " number of subforce_eval specified.",&
            usage="group_partition  2 2 4 2 4 ", type_of_var=integer_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NGROUPS",variants=(/"NGROUP"/),&
            description="Gives the wanted number of groups. If not specified the number"//&
            " of groups is set to the number of subforce_eval defined.",&
            usage="ngroups 4", type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       ! Double force_eval
       CALL section_create(subsection,name="LINEAR",&
            description="Linear combination between two force_eval:  F= lambda F1 + (1-lambda) F2",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       CALL keyword_create(keyword, name="LAMBDA",&
            description="Specify the mixing parameter lambda in the formula:",&
            usage="lambda <REAL>", type_of_var=real_t,&
            required=.TRUE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="COUPLING",&
            description="Coupling between two force_eval: E=(E1+E2 - sqrt((E1-E2)**2+4*H12**2))/2",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       CALL keyword_create(keyword, name="COUPLING_PARAMETER",&
            description="Coupling parameter H12 used in the coupling",&
            usage="COUPLING_PARAMETER <REAL>", type_of_var=real_t,&
            required=.TRUE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="RESTRAINT",&
            description="Restraint between two force_eval: E = E1 + k*(E1-E2-t)**2",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       CALL keyword_create(keyword, name="RESTRAINT_TARGET",&
            description="Target value of the restraint (t) ",&
            usage="RESTRAINT_TARGET <REAL>", type_of_var=real_t,&
            required=.TRUE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTRAINT_STRENGTH",&
            description="Strength of the restraint (k) in "//&
            "k*(E1-E2-t)**2" ,&
            usage="RESTRAINT_STRENGTH <REAL>", type_of_var=real_t,&
            required=.TRUE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)
       
       ! Multiple force_eval
       CALL section_create(subsection,name="GENERIC",&
            description="User driven coupling between two or more force_eval.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       CALL keyword_create(keyword, name="MIXING_FUNCTION",&
            description="Specifies the mixing functional form in mathematical notation.",&
            usage="MIXING_FUNCTION (E1+E2-LOG(E1/E2))", required=.TRUE., type_of_var=char_t,&
            n_var=-1, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="VARIABLES",&
            description="Defines the variables of the functional form. To allow an efficient"//&
            " mapping the order of the energy variables will be considered identical to the"//&
            " order of the force_eval in the force_eval_order list.",&
            usage="VARIABLES x", required=.TRUE., type_of_var=char_t,&
            n_var=-1, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PARAMETERS",&
            description="Defines the parameters of the functional form",&
            usage="PARAMETERS a b D", required=.TRUE., type_of_var=char_t,&
            n_var=-1, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="VALUES",&
            description="Defines the values of  parameter of the functional form",&
            usage="VALUES ", required=.TRUE., type_of_var=real_t,&
            n_var=-1, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DX",&
            description="Parameter used for computing the derivative with the Ridders method.",&
            usage="DX <REAL>", default_r_val=0.1_dp, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ERROR_LIMIT",&
            description="Checks that the error in computing the derivative is not larger than "//&
            "the value set. In case prints a warning message.",&
            usage="ERROR_LIMIT <REAL>", default_r_val=1.0E-12_dp, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       ! Mapping of atoms
       NULLIFY(sub2section, sub3section)
       CALL section_create(subsection,name="MAPPING",&
            description="Defines the mapping of atoms for the different force_eval with the mixed force_eval."//&
            " The default is to have a mapping 1-1 between atom index (i.e. all force_eval share the same"//&
            " geometrical structure). The mapping is based on defining fragments and the mapping the "//&
            " fragments between the several force_eval and the mixed force_eval",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       ! Mixed force_eval
       CALL section_create(sub2section,name="FORCE_EVAL_MIXED",&
            description="Defines the fragments for the mixed force_eval (reference)",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL section_create(sub3section,name="FRAGMENT",&
            description="Fragment definition",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="Defines the index of the fragment defined",&
            usage="<INTEGER>",  type_of_var=integer_t, n_var=1, required=.TRUE.,error=error)
       CALL section_add_keyword(sub3section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Starting and ending atomic index defining one fragment must be provided",& 
            usage="<INTEGER> <INTEGER>", type_of_var=integer_t, n_var=2, repeats=.TRUE.,&
            required=.TRUE.,error=error)
       CALL section_add_keyword(sub3section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(sub2section,sub3section,error=error)
       CALL section_release(sub3section,error=error)
       CALL section_add_subsection(subsection,sub2section,error=error)
       CALL section_release(sub2section,error=error)

       ! All other force_eval
       CALL section_create(sub2section,name="FORCE_EVAL",&
            description="Defines the fragments and the mapping for each force_eval (an integer index (ID) "//&
            "needs to be provided as parameter)",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="DEFINE_FRAGMENTS",&
            description="Specify the fragments definition of the force_eval through the fragments of the"//&
            " force_eval_mixed. This avoids the pedantic definition of the fragments for the force_eval,"//&
            " assuming the order of the fragments for the specified force_eval is the same as the sequence "//&
            " of integers provided. Easier to USE should be preferred to the specification of the single fragments.",&
            usage="DEFINE_FRAGMENTS <INTEGER> .. <INTEGER>", type_of_var=integer_t, n_var=-1,&
            error=error)
       CALL section_add_keyword(sub2section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="Defines the index of the force_eval for which fragments and mappings are provided",&
            usage="<INTEGER>",  type_of_var=integer_t, n_var=1, required=.TRUE.,error=error)
       CALL section_add_keyword(sub2section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_create(sub3section,name="FRAGMENT",&
            description="Fragment definition",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="Defines the index of the fragment defined",&
            usage="<INTEGER>",  type_of_var=integer_t, n_var=1, required=.FALSE.,error=error)
       CALL section_add_keyword(sub3section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Starting and ending atomic index defining one fragment must be provided",& 
            usage="<INTEGER> <INTEGER>", type_of_var=integer_t, n_var=2, repeats=.FALSE.,&
            required=.TRUE.,error=error)
       CALL section_add_keyword(sub3section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAP",&
            description="Provides the index of the fragment of the MIXED force_eval mapped on the"//&
            " locally defined fragment.",& 
            usage="MAP <INTEGER>", type_of_var=integer_t, n_var=1, repeats=.FALSE.,&
            required=.TRUE.,error=error)
       CALL section_add_keyword(sub3section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(sub2section,sub3section,error=error)
       CALL section_release(sub3section,error=error)
       CALL section_add_subsection(subsection,sub2section,error=error)
       CALL section_release(sub2section,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="controls the printing of ewald setup",&
            print_level=low_print_level,each=(/1/),add_last=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)
       
    END IF

  END SUBROUTINE create_mix_section


END MODULE input_cp2k
