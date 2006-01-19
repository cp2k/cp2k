!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
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
  USE input_cp2k_motion, only: create_motion_section
  USE input_cp2k_subsys, only: create_subsys_section
  USE input_cp2k_mm, only: create_mm_section
  USE input_cp2k_qmmm, only: create_qmmm_section
  USE input_cp2k_dft, ONLY: create_dft_section, create_bsse_section, &
       create_properties_section
  USE input_cp2k_tbmd, ONLY: create_tbmd_section

#include "cp_common_uses.h"
  IMPLICIT NONE
  PRIVATE
  
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='input_cp2k'
  
  PUBLIC :: create_cp2k_root_section, check_cp2k_input, parsed_cp2k_input,&
       create_global_section
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

       CALL create_farming_section(section,error=error)
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
    TYPE(section_type), pointer :: print_key

    failure=.false.
    nullify(print_key)
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
         description="Which program should be run to compute forces",&
         usage="PROGRAM_NAME <STRING>",&
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
         usage="PROJECT_NAME <STRING>",&
         default_c_val="PROJECT",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="OUTPUT_FILE_NAME",&
         description="Name of the output file. "//&
         "Relevant only if automatically started (through farming for example). "//&
         "If empty uses the project name as basis for it.",&
         usage="OUTPUT_FILE_NAME {filename}",default_c_val="", error=error)
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

    CALL keyword_create(keyword, name="SEED",&
         description="Seed of the random generator used for the velocities. "//&
         "Should be negative",&
         usage="SEED -2", default_i_val=-1, error=error)
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

    CALL cp_print_key_section_create(print_key, "print", description=&
         "controls the printing of physical and mathematical constants", &
         print_level=high_print_level,filename="__STD_OUT__",&
         error=error)
    CALL keyword_create(keyword,"physcon",&
         description="if the printkey is active prints the physical constants",&
         default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword,"spherical_harmonics",&
         description="if the printkey is active prints the spherical harmonics",&
         default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword,"rgn_matrices",&
         description="if the printkey is active prints the random number generator matrices",&
         default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key, "performance", description=&
         "controls the printing of physical and mathematical constants", &
         print_level=high_print_level,filename="__STD_OUT__",&
         error=error)
    CALL keyword_create(keyword,"HPM",&
         description="if the printkey is active controls the printing of the hpm info (if linked as"//&
         " EXTERNAL library)",&
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

    CALL cp_print_key_section_create(print_key, "parallel", description=&
         "controls the printing of physical and mathematical constants", &
         print_level=high_print_level,filename="__STD_OUT__",&
         error=error)
    CALL keyword_create(keyword,"mpi",&
         description="if the printkey is active controls the printing of"//&
         " the hpm info (IF linked as EXTERNAL library)",&
         default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword,"blacs_info",&
         description="if the printkey is active controls the printing of"//&
         " the blacs info (IF linked as EXTERNAL library)",&
         default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
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
!***************************************************************************

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
    TYPE(section_type), pointer :: section
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_test_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword

    call section_create(section,name="test",&
         description="tests to perform",&
         n_keywords=7, n_subsections=0, repeats=.FALSE., required=.true.,&
         error=error)

    nullify(keyword)
    CALL keyword_create(keyword, name="memory",&
         description="tries to allocate the given amount of memory ",&
         usage="BLACS_GRID SQUARE",default_r_val=256.e6_dp, error=error)
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

    CALL keyword_create(keyword, name="random_number_generator",variants=(/"rgn"/),&
         description=" Tests the parallel random number generator",&
         usage="random_number_generator",default_l_val=.FALSE.,&
         lone_keyword_l_val=.true.,error=error)
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

  END SUBROUTINE create_test_section
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
         n_keywords=0, n_subsections=8, repeats=.FALSE., required=.TRUE.,&
         error=error)


    nullify(subsection)
    CALL create_dft_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_tbmd_section(subsection, error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection, error=error)

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

    CALL create_f_env_print_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

  END SUBROUTINE create_force_eval_section
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
  TYPE(section_type), pointer :: section
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='create_farming_section',&
        routineP=moduleN//':'//routineN
  TYPE(section_type), pointer :: sub_section
  TYPE(keyword_type), pointer :: keyword

  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  CALL section_create(section,name="farming",&
         description="describes a farming job, in which multiple inputs are executed",&
         repeats=.FALSE., required=.TRUE., error=error)
  nullify(keyword)

  CALL keyword_create(keyword, name="group_size",&
         description="gives the wanted size of a working group.",&
         usage="group_size 2", default_i_val=8, error=error)
  CALL section_add_keyword(section,keyword,error=error)
  CALL keyword_release(keyword,error=error)

  CALL keyword_create(keyword, name="ngroups",variants=(/"ngroup"/),&
         description="gives the wanted number of groups. Overrides GROUP_SIZE.",&
         usage="ngroups 4", type_of_var=integer_t, error=error)
  CALL section_add_keyword(section,keyword,error=error)
  CALL keyword_release(keyword,error=error)

  CALL keyword_create(keyword, name="group_partition",&
         description="gives the exact number of processors for each group, ovverides both NGROUPS and GROUP_PARTITION",&
         usage="ngroups 4", type_of_var=integer_t, n_var=-1, error=error)
  CALL section_add_keyword(section,keyword,error=error)
  CALL keyword_release(keyword,error=error)

  CALL keyword_create(keyword, name="restart",&
         description="if you are restarting a farming job (and should pick up where the previous left off)",&
         usage="RESTART",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,&
         error=error)
  CALL section_add_keyword(section,keyword,error=error)
  CALL keyword_release(keyword,error=error)

  CALL keyword_create(keyword, name="max_step",variants=(/"max_steps"/),&
         description="maximum number of steps done by a processor",&
         usage="max_step 4", default_i_val=HUGE(0)/32768, error=error)
  CALL section_add_keyword(section,keyword,error=error)
  CALL keyword_release(keyword,error=error)

  CALL keyword_create(keyword, name="cycle",&
         description="If farming should cycle between the inputs" ,&
         usage="CYCLE",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,&
         error=error)
  CALL section_add_keyword(section,keyword,error=error)
  CALL keyword_release(keyword,error=error)

  nullify(sub_section)
  CALL section_create(sub_section,name="job",&
         description="description of the jobs to be executed",&
         repeats=.true., required=.TRUE., error=error)

  CALL keyword_create(keyword,name="directory",&
       description="the directory in which the job should be executed",&
       usage="directory /my/path",&
       required=.TRUE.,default_c_val=".",error=error)
  CALL section_add_keyword(sub_section,keyword,error=error)
  CALL keyword_release(keyword,error=error)

  CALL keyword_create(keyword,name="input_file_name",&
       description="the filename of the input file",&
       usage="input_file my_input.inp",&
       required=.TRUE.,default_c_val="input.inp",error=error)
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
  TYPE(section_type), pointer :: section
  type(cp_error_type), optional, intent(inout) :: error

  logical :: failure
  character(len=*), parameter :: routineN='create_f_env_print_section',&
       routineP=moduleN//':'//routineN
  TYPE(section_type), pointer :: print_key
  TYPE(keyword_type), pointer :: keyword

  failure=.FALSE.

  NULLIFY(print_key,keyword)
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL section_create(section,name="print",&
          description="propreties that you want to output and that are common to all methods",&
          n_keywords=0, n_subsections=5, repeats=.FALSE., required=.TRUE.,&
          error=error)

     CALL cp_print_key_section_create(print_key,"trajectory",&
          description="controls the output of the trajectory",&
          print_level=low_print_level, filename="",error=error)
     CALL keyword_create(keyword, name="format",&
          description="Specifies the format of the output file for the trajectory.",&
          usage="FORMAT (ATOMIC|XMOL|DCD)",&
          default_i_val=dump_xmol,&
          enum_c_vals=output_name_types,&
          enum_i_vals=output_name_index,&
          error=error)
     CALL section_add_keyword(print_key,keyword,error=error)
     CALL keyword_release(keyword,error=error)    
     CALL section_add_subsection(section,print_key,error=error)
     CALL section_release(print_key,error=error)

     CALL cp_print_key_section_create(print_key,"velocities",&
          description="controls the output of the velocities",&
          print_level=low_print_level, filename="",error=error)
     CALL keyword_create(keyword, name="format",&
          description="Specifies the format of the output file for the trajectory.",&
          usage="FORMAT (ATOMIC|XMOL|DCD)",&
          default_i_val=dump_xmol,&
          enum_c_vals=output_name_types,&
          enum_i_vals=output_name_index,&
          error=error)
     CALL section_add_keyword(print_key,keyword,error=error)
     CALL keyword_release(keyword,error=error)    
     CALL section_add_subsection(section,print_key,error=error)
     CALL section_release(print_key,error=error)

     CALL cp_print_key_section_create(print_key,"forces",&
          description="controls the output of the forces",&
          print_level=high_print_level, filename="",error=error)
     CALL keyword_create(keyword, name="format",&
          description="Specifies the format of the output file for the trajectory.",&
          usage="FORMAT (ATOMIC|XMOL|DCD)",&
          default_i_val=dump_xmol,&
          enum_c_vals=output_name_types,&
          enum_i_vals=output_name_index,&
          error=error)
     CALL section_add_keyword(print_key,keyword,error=error)
     CALL keyword_release(keyword,error=error)    
     CALL section_add_subsection(section,print_key,error=error)
     CALL section_release(print_key,error=error)

     CALL cp_print_key_section_create(print_key,"energy",&
          description="controls the output the ener file",&
          print_level=debug_print_level, filename="",error=error)
     CALL section_add_subsection(section,print_key,error=error)
     CALL section_release(print_key,error=error)

     CALL cp_print_key_section_create(print_key,"atomic_coordinates",&
          description="controls the output of the atomic coordinates (check if trajectory is not what you want)",&
          print_level=low_print_level, filename="__STD_OUT__",error=error)
     CALL section_add_subsection(section,print_key,error=error)
     CALL section_release(print_key,error=error)

     CALL cp_print_key_section_create(print_key,"interatomic_distances",&
          description="controls the output of the interatomic distances",&
          print_level=debug_print_level, filename="__STD_OUT__",error=error)
     CALL section_add_subsection(section,print_key,error=error)
     CALL section_release(print_key,error=error)

     CALL cp_print_key_section_create(print_key,"cell",&
          description="controls the output of the cell parameters",&
          print_level=medium_print_level, filename="",error=error)
     CALL section_add_subsection(section,print_key,error=error)
     CALL section_release(print_key,error=error)

     CALL cp_print_key_section_create(print_key,"kinds",&
          description="controls the output of information on the kinds",&
          print_level=medium_print_level, filename="",error=error)
     CALL keyword_create(keyword, name="kind_information",&
          description="If the printkey is activated controls the printing of generic kind information",&
          default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
     CALL section_add_keyword(print_key,keyword,error=error)
     CALL keyword_release(keyword,error=error)    
     CALL keyword_create(keyword, name="gth_potential",&
          description="If the printkey is activated controls the printing of the gth_potential information",&
          default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
     CALL section_add_keyword(print_key,keyword,error=error)
     CALL keyword_release(keyword,error=error)    
     CALL keyword_create(keyword, name="basis_set",&
          description="If the printkey is activated controls the printing of basis set information",&
          default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
     CALL section_add_keyword(print_key,keyword,error=error)
     CALL keyword_release(keyword,error=error)    
     CALL section_add_subsection(section,print_key,error=error)
     CALL section_release(print_key,error=error)

     CALL cp_print_key_section_create(print_key,"radii",&
          description="controls the output of radii information",&
          print_level=high_print_level, filename="",error=error)
     CALL keyword_create(keyword, name="core_charges_radii",&
          description="If the printkey is activated controls the printing of the radii of the core charges",&
          default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
     CALL section_add_keyword(print_key,keyword,error=error)
     CALL keyword_release(keyword,error=error)    
     CALL keyword_create(keyword, name="pgf_radii",&
          description="If the printkey is activated controls the printing of the core gaussian radii",&
          default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
     CALL section_add_keyword(print_key,keyword,error=error)
     CALL keyword_release(keyword,error=error)    
     CALL keyword_create(keyword, name="set_radii",&
          description="If the printkey is activated controls the printing of the set_radii",&
          default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
     CALL section_add_keyword(print_key,keyword,error=error)
     CALL keyword_release(keyword,error=error)    
     CALL keyword_create(keyword, name="gth_ppl_radii",&
          description="If the printkey is activated controls the printing of the "//&
          "gth pseudo potential local radii",&
          default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
     CALL section_add_keyword(print_key,keyword,error=error)
     CALL keyword_release(keyword,error=error)    
     CALL keyword_create(keyword, name="gth_ppnl_radii",&
          description="If the printkey is activated controls the printing of the "//&
          "gth pseudo potential non local radii",&
          default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
     CALL section_add_keyword(print_key,keyword,error=error)
     CALL keyword_release(keyword,error=error)    
     CALL keyword_create(keyword, name="gapw_prj_radii",&
          description="If the printkey is activated controls the printing of the gapw projector radii",&
          default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
     CALL section_add_keyword(print_key,keyword,error=error)
     CALL keyword_release(keyword,error=error)    
     CALL section_add_subsection(section,print_key,error=error)
     CALL section_release(print_key,error=error)

     CALL cp_print_key_section_create(print_key,"total_numbers",&
          description="controls the printing of the total number of atoms, kinds,...",&
          print_level=medium_print_level, filename="__STD_OUT__",error=error)
     CALL section_add_subsection(section,print_key,error=error)
     CALL section_release(print_key,error=error)

     CALL cp_print_key_section_create(print_key,"distribution",&
          description="controls the printing of the distribution of molecules, atoms, matrix blocks,...",&
          print_level=medium_print_level, filename="__STD_OUT__",error=error)
     CALL section_add_subsection(section,print_key,error=error)
     CALL section_release(print_key,error=error)

  END IF
end subroutine create_f_env_print_section
!***************************************************************************

END MODULE input_cp2k
