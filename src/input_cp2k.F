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
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE input_constants
  USE input_cp2k_dft,                  ONLY: create_bsse_section,&
                                             create_dft_section,&
                                             create_properties_section
  USE input_cp2k_eip,                  ONLY: create_eip_section
  USE input_cp2k_mm,                   ONLY: create_mm_section
  USE input_cp2k_motion,               ONLY: create_motion_section
  USE input_cp2k_qmmm,                 ONLY: create_qmmm_section
  USE input_cp2k_subsys,               ONLY: create_subsys_section
  USE input_cp2k_tbmd,                 ONLY: create_tbmd_section
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type,&
                                             section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type
  USE input_val_types,                 ONLY: char_t,&
                                             integer_t,&
                                             no_t
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
       create_global_section
!!***
!****************************************************************************
CONTAINS

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_cp2k_root_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section

    CALL timeset(routineN,handle)
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

       CALL create_motion_section(section,error=error)
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

       NULLIFY(keyword)
       CALL keyword_create(keyword,name="_DEFAULT_KEYWORD_",&
            description="ignores lines outside sections",&
            usage=" ", type_of_var=no_t,error=error)
       CALL section_add_keyword(root_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
    CALL timestop(handle)

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_global_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.
    NULLIFY(print_key)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="global",&
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

    CALL keyword_create(keyword, name="FFT_LIBRARY",&
         variants=s2a("FFT_LIB","FFTLIB"),&
         description="which fft library to use to perform the ffts",&
         usage="FFT_LIBRARY FFTW",&
         default_c_val="FFTSG",error=error)
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
         print_level=low_print_level,filename="__STD_OUT__",&
         error=error)
    CALL keyword_create(keyword,"BASIC_DATA_TYPES",&
         description="Controls the printing of the basic data types.",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword,"physcon",&
         description="if the printkey is active prints the physical constants",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword,"spherical_harmonics",&
         description="if the printkey is active prints the spherical harmonics",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword,"rng_matrices",&
         description="if the printkey is active prints the random number generator matrices",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_test_section', &
      routineP = moduleN//':'//routineN

    TYPE(keyword_type), POINTER              :: keyword

    CALL section_create(section,name="test",&
         description="tests to perform",&
         n_keywords=7, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
         error=error)

    NULLIFY(keyword)
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
         description=" Tests the parallel random number generator",&
         usage="random_number_generator",default_l_val=.FALSE.,&
         lone_keyword_l_val=.TRUE.,error=error)
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
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_force_eval_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="force_eval",&
         description="parameters needed to calculate energy and forces and"//&
         " describe the system you want to analyze.",&
         n_keywords=0, n_subsections=9, repeats=.FALSE., required=.TRUE.,&
         error=error)

    NULLIFY(subsection)

    CALL create_dft_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_tbmd_section(subsection, error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection, error=error)

    CALL create_eip_section(subsection, error)
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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

       CALL keyword_create(keyword, name="EXTERNAL_FILE",&
            description="Specify tha external input file to read",&
            default_lc_val=" ", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART_COUNTERS",&
            description="Restarts the counters in MD schemes",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_POS",&
            description="Takes the positions from the external file",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_VEL",&
            description="Takes the velocities from the external file",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
      CALL keyword_create(keyword, name="RESTART_RANDOMG",&
            description="Restarts the random number generator from the external file",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_KG_COEFF_POS",&
            description="Restarts the positions of coefficents of kg from the external file",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_KG_COEFF_VEL",&
            description="Restarts the velocity of coefficents of kg from the external file",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_KG_NOSE",&
            description="Restarts the nose thermostats of the kg coefficents "//&
            "from the EXTERNAL file",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_BAROSTAT",&
            description="Restarts the barostat from the external file",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_BAROSTAT_NOSE",&
            description="Restarts the barostat thermostat from the external file",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_NOSE",&
            description="Restarts the nose thermostats of the particles "//&
            "from the EXTERNAL file",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_CELL",&
            description="Restarts the cell (and cell_ref) "//&
            "from the EXTERNAL file",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_METADYNAMICS",&
            description="Restarts hills from a previous metadynamics run "//&
            "from the EXTERNAL file",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE.,error=error)
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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

  NULLIFY(sub_section)
  CALL section_create(sub_section,name="job",&
         description="description of the jobs to be executed",&
         repeats=.TRUE., required=.TRUE., error=error)

  CALL keyword_create(keyword,name="directory",&
       description="the directory in which the job should be executed",&
       usage="directory /my/path",&
       required=.TRUE.,default_lc_val=".",error=error)
  CALL section_add_keyword(sub_section,keyword,error=error)
  CALL keyword_release(keyword,error=error)

  CALL keyword_create(keyword,name="input_file_name",&
       description="the filename of the input file",&
       usage="input_file my_input.inp",&
       required=.TRUE.,default_lc_val="input.inp",error=error)
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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

     CALL cp_print_key_section_create(print_key,"PW_GRID_INFORMATION",&
          description="Controls the printing of information regarding the pw_grid structure.",&
          print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
          error=error)
     CALL section_add_subsection(section,print_key,error=error)
     CALL section_release(print_key,error=error)

     CALL cp_print_key_section_create(print_key,"RS_GRID_INFORMATION",&
          description="Controls the printing of information regarding the rs_grid structure.",&
          print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
          error=error)
     CALL section_add_subsection(section,print_key,error=error)
     CALL section_release(print_key,error=error)

     CALL cp_print_key_section_create(print_key,"total_numbers",&
          description="controls the printing of the total number of atoms, kinds,...",&
          print_level=medium_print_level, filename="__STD_OUT__",error=error)
     CALL section_add_subsection(section,print_key,error=error)
     CALL section_release(print_key,error=error)

     CALL cp_print_key_section_create(print_key,"distribution",&
          description="controls the printing of the distribution of molecules, atoms, ...",&
          print_level=medium_print_level, filename="__STD_OUT__",error=error)
     CALL section_add_subsection(section,print_key,error=error)
     CALL section_release(print_key,error=error)

     CALL cp_print_key_section_create(print_key,"distribution2d",&
          description="controls the printing of the distribution of matrix blocks,...",&
          print_level=high_print_level, filename="__STD_OUT__",error=error)
     CALL section_add_subsection(section,print_key,error=error)
     CALL section_release(print_key,error=error)

     CALL cp_print_key_section_create(print_key,"distribution1d",&
          description="Each node prints out its distribution info...",&
          print_level=high_print_level, filename="__STD_OUT__",error=error)
     CALL section_add_subsection(section,print_key,error=error)
     CALL section_release(print_key,error=error)

  END IF
END SUBROUTINE create_f_env_print_section
!***************************************************************************

END MODULE input_cp2k
