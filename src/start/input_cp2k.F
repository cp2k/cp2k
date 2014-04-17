!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief builds the input structure for cp2k
!> \par History
!>      06.2004 created [fawzi]
!> \author fawzi
! *****************************************************************************
MODULE input_cp2k
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE dbcsr_api,                       ONLY: dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8
  USE f77_blas
  USE input_constants
  USE input_cp2k_atom,                 ONLY: create_atom_section
  USE input_cp2k_force_eval,           ONLY: create_force_eval_section
  USE input_cp2k_global,               ONLY: create_global_section
  USE input_cp2k_motion,               ONLY: create_motion_section
  USE input_cp2k_rsgrid,               ONLY: create_rsgrid_section
  USE input_cp2k_vib,                  ONLY: create_vib_section
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_optimize_basis,            ONLY: create_optimize_basis_section
  USE input_optimize_input,            ONLY: create_optimize_input_section
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: char_t,&
                                             integer_t,&
                                             lchar_t,&
                                             logical_t,&
                                             real_t
  USE kinds,                           ONLY: dp
  USE pw_grids,                        ONLY: do_pw_grid_blocked_false,&
                                             do_pw_grid_blocked_free,&
                                             do_pw_grid_blocked_true
  USE string_utilities,                ONLY: newline,&
                                             s2a
  USE swarm_input,                     ONLY: create_swarm_section
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k'

  PUBLIC :: create_cp2k_root_section

CONTAINS

! *****************************************************************************
!> \brief creates the input structure of the file used by cp2k
!> \param root_section the input structure to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE create_cp2k_root_section(root_section,error)
    TYPE(section_type), POINTER              :: root_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_cp2k_root_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: section

    CALL timeset(routineN,handle)
    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(root_section),cp_failure_level,routineP,error,failure)
    IF (.NOT.  failure) THEN
       CALL section_create(root_section,name="__ROOT__",&
            description="input file of cp2k",n_keywords=0, n_subsections=10,&
            repeats=.FALSE., required=.FALSE., supported_feature=.TRUE., error=error)
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

       CALL create_optimize_input_section(section,error=error)
       CALL section_add_subsection(root_section,section,error=error)
       CALL section_release(section,error=error)

       CALL create_optimize_basis_section(section,error=error)
       CALL section_add_subsection(root_section,section,error=error)
       CALL section_release(section,error=error)

       CALL create_swarm_section(section,error=error)
       CALL section_add_subsection(root_section,section,error=error)
       CALL section_release(section,error=error)

       CALL create_ext_restart_section(section,error=error)
       CALL section_add_subsection(root_section,section,error=error)
       CALL section_release(section,error=error)

       CALL create_vib_section(section,error=error)
       CALL section_add_subsection(root_section,section,error=error)
       CALL section_release(section,error=error)

       CALL create_atom_section(section,error=error)
       CALL section_add_subsection(root_section,section,error=error)
       CALL section_release(section,error=error)
    END IF
    CALL timestop(handle)

  END SUBROUTINE create_cp2k_root_section


! *****************************************************************************
!> \brief section with the tests of the libraries or external code that cp2k uses
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE create_test_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_test_section', &
      routineP = moduleN//':'//routineN

    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    CALL section_create(section,name="TEST",&
         description="Tests to perform on the supported libraries.",&
         n_keywords=7, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
         error=error)

    NULLIFY(keyword, print_key)
    CALL keyword_create(keyword, name="MEMORY",&
         description="Set the maximum amount of memory allocated for a given test (in bytes)",&
         usage="MEMORY <REAL>",default_r_val=256.e6_dp, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="COPY",&
         description="Tests the performance to copy two vectors."//&
         "The results of these tests allow to determine the size of the cache "//&
         "of the CPU. This can be used to optimize the performance of the"//&
         "FFTSG library. Tests are repeated the given number of times.",&
         usage="copy 10",default_i_val=0,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MATMUL",&
         description="Tests the performance of different kinds of matrix matrix "//&
         "multiply kernels for the F95 INTRINSIC matmul. Matrices up to 2**N+1 will be tested. ",&
         usage="matmul 10",default_i_val=0,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="DGEMM",&
         description="Tests the performance of different kinds of matrix matrix "//&
         "multiply kernels for the BLAS INTRINSIC DGEMM. Matrices up to 2**N+1 will be tested. ",&
         usage="DGEMM 10",default_i_val=0,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="FFT",&
         description="Tests the performance of all available FFT libraries for "//&
         "3D FFTs Tests are repeated the given number of times.",&
         usage="fft 10",default_i_val=0,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="ERI",&
         description="Tests the performance and correctness of ERI libraries ",&
         usage="eri 1",default_i_val=0,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="CLEBSCH_GORDON",variants=(/"CLEBSCH"/),&
         description="Tests the Clebsch-Gordon Coefficients. "//&
         "Tests are repeated the given number of times.",&
         usage="clebsch_gordon 10",default_i_val=0,error=error)

    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MPI",&
         description="Tests mpi, quickly adapted benchmark code,"//&
         "will ONLY work on an even number of CPUs. comm is the relevant, "//&
         "initialized communicator. This test will produce messages "//&
         "of the size 8*10**requested_size, where requested_size is the value "//&
         "given to this keyword",&
         usage="mpi 6",default_i_val=0,error=error)

    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="RANDOM_NUMBER_GENERATOR",variants=(/"rng"/),&
         description=" Tests the parallel random number generator (RNG)",&
         usage="rng 1000000",default_i_val=0,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL cp_print_key_section_create(print_key,"GRID_INFORMATION",&
         description="Controls the printing of information regarding the PW and RS grid structures"//&
         " (ONLY for TEST run).",&
         print_level=medium_print_level,filename="__STD_OUT__",&
         error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
         description="controls the printing of tests output",&
         print_level=silent_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
         error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    NULLIFY(subsection)
    CALL create_rs_pw_transfer_section(subsection,error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

    CALL create_eigensolver_section(subsection,error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

    CALL create_pw_transfer_section(subsection,error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

    CALL create_cp_fm_gemm_section(subsection,error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

    CALL create_cp_dbcsr_section(subsection,error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

  END SUBROUTINE create_test_section


! *****************************************************************************
!> \brief section to setup debugging parameter
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_debug_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_debug_section', &
      routineP = moduleN//':'//routineN

    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    CALL section_create(section,name="DEBUG",&
         description="Section to setup parameters for debug runs.",&
         n_keywords=7, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
         error=error)

    NULLIFY (keyword,print_key)

    CALL keyword_create(keyword, name="DEBUG_FORCES",&
         description="Activates the debugging of the atomic forces",&
         usage="DEBUG_FORCES {LOGICAL}",default_l_val=.TRUE.,&
         lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="DEBUG_STRESS_TENSOR",&
         description="Activates the debugging of the stress tensor",&
         usage="DEBUG_STRESS_TENSOR {LOGICAL}",default_l_val=.TRUE.,&
         lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="DX",&
         description="Increment for the calculation of the numerical derivatives",&
         usage="DX {REAL}",default_r_val=0.001_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="EPS_NO_ERROR_CHECK",&
         description="The mismatch between the numerical and the "//&
                     "analytical value is not checked for analytical "//&
                     "values smaller than this threshold value",&
         usage="EPS_NO_ERROR_CHECK {REAL}",default_r_val=1.0E-5_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="STOP_ON_MISMATCH",&
         description="Stop the debug run when a mismatch between the numerical and "//&
                     "the analytical value is detected",&
         usage="STOP_ON_MISMATCH {LOGICAL}",default_l_val=.TRUE.,&
         lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
         description="Controls the printing of the DEBUG specific output",&
         print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
         error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

  END SUBROUTINE create_debug_section


! *****************************************************************************
!> \brief creates the multiple force_eval section
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
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
         supported_feature=.TRUE.,error=error)

    NULLIFY(keyword)
    CALL keyword_create(keyword, name="FORCE_EVAL_ORDER",&
         description='Specify the orders of the different force_eval. When using a MIXED force_eval'//&
         " this does not need to be specified in this list, because it that takes into account only the real"//&
         " energy contributions",&
         usage="FORCE_EVAL_ORDER <INTEGER> .. <INTEGER>", type_of_var=integer_t, n_var=-1,&
         default_i_vals=(/1/), supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MULTIPLE_SUBSYS",&
         description="Specify if force_eval have different subsys. In case they share the same subsys,"//&
         " it needs to be specified only in the MIXED force_eval (if using MIXED) or"//&
         " in the force_eval corresponding to first force_eval of the previous order (when not using MIXED).",&
         default_l_val=.FALSE., lone_keyword_l_val=.TRUE., supported_feature=.TRUE., error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

   END SUBROUTINE create_multi_force_section


! *****************************************************************************
!> \brief Creates the exteranal restart section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
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
            supported_feature=.TRUE.,error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="RESTART_FILE_NAME",variants=(/"EXTERNAL_FILE"/),&
            description="Specifies the name of restart file (or any other input file)"//&
                        " to be read. Only fields relevant to a restart will be used"//&
                        " (unless switched off with the keywords in this section)",&
            default_lc_val=" ", supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BINARY_RESTART_FILE_NAME",&
                           variants=(/"BINARY_RESTART_FILE"/),&
                           description="Specifies the name of an additional restart file "//&
                           "from which selected input sections are read in binary format "//&
                           "(see SPLIT_RESTART_FILE).",&
                           default_lc_val="",&
                           supported_feature=.TRUE.,&
                           error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART_DEFAULT",&
            description="This keyword controls the default value for all possible "//&
            " restartable keywords, unless explicitly defined. For example setting"//&
            " this keyword to .FALSE. does not restart any quantity. If, at the "//&
            " same time, one keyword is set to .TRUE. only that quantity will be"//&
            " restarted.", default_l_val=.TRUE., supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_COUNTERS",&
            description="Restarts the counters in MD schemes",&
            type_of_var=logical_t, lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_POS",&
            description="Takes the positions from the external file",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_VEL",&
            description="Takes the velocities from the external file",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_RANDOMG",&
            description="Restarts the random number generator from the external file",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART_SHELL_POS",&
            description="Takes the positions of the shells from the external file (only if shell-model)",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_CORE_POS",&
            description="Takes the positions of the cores from the external file (only if shell-model)",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_OPTIMIZE_INPUT_VARIABLES",&
            description="Restart with the optimize input variables",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART_SHELL_VELOCITY",&
            description="Takes the velocities of the shells from the external file (only if shell-model)",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_CORE_VELOCITY",&
            description="Takes the velocities of the shells from the external file (only if shell-model)",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_BAROSTAT",&
            description="Restarts the barostat from the external file",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_BAROSTAT_THERMOSTAT",&
            description="Restarts the barostat thermostat from the external file",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_SHELL_THERMOSTAT",&
            description="Restarts the shell thermostat from the external file",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_THERMOSTAT",&
            description="Restarts the nose thermostats of the particles "//&
            "from the EXTERNAL file",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_CELL",&
            description="Restarts the cell (and cell_ref) "//&
            "from the EXTERNAL file",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_METADYNAMICS",&
            description="Restarts hills from a previous metadynamics run "//&
            "from the EXTERNAL file",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_WALKERS",&
            description="Restarts walkers informations from a previous metadynamics run "//&
            "from the EXTERNAL file",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_BAND",&
            description="Restarts positions and velocities of the Band.",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_QMMM",&
            description="Restarts the following specific QMMM info: translation vectors.",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_CONSTRAINT",&
            description="Restarts constraint section. It's necessary when doing restraint"//&
            " calculation to have a perfect energy conservation. For constraints only it's"//&
            " use is optional.",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_BSSE",&
            description="Restarts information for BSSE calculations.",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_DIMER",&
            description="Restarts information for DIMER geometry optimizations.",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_AVERAGES",&
            description="Restarts information for AVERAGES.",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_RTP",&
            description="Restarts information for REAL TIME PROPAGATION and EHRENFEST DYNAMICS.",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="CUSTOM_PATH",&
            description="Restarts the given path from the EXTERNAL file. Allows a major flexibility for restarts.",&
            type_of_var=char_t,repeats=.TRUE., supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! PIMD
       CALL keyword_create(keyword, name="RESTART_PINT_POS",&
            description="Restart bead positions from PINT%BEADS%COORD.",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_PINT_VEL",&
            description="Restart bead velocities from PINT%BEADS%VELOCITY.",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_PINT_NOSE",&
            description="Restart Nose thermostat for beads from PINT%NOSE.",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_PINT_GLE",&
            description="Restart GLE thermostat for beads from PINT%GLE.",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       ! PIMC
       CALL keyword_create(keyword, name="RESTART_HELIUM_POS",&
            description="Restart helium positions from PINT%HELIUM%COORD.",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_HELIUM_PERMUTATION",&
            description="Restart helium permutation state from PINT%HELIUM%PERM.",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_HELIUM_FORCE",&
            description="Restart helium forces exerted on the solute from PINT%HELIUM%FORCE.",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_HELIUM_RNG",&
            description="Restarts helium random number generators from PINT%HELIUM%RNG_STATE.",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RESTART_HELIUM_DENSITIES",&
            description="Restarts helium density distributions from PINT%HELIUM%RHO.",&
            type_of_var=logical_t,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,&
            default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_ext_restart_section


! *****************************************************************************
!> \brief creates the farming section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE create_farming_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_farming_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, sub_section

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="farming",&
         description="Describes a farming job, in which multiple inputs are executed."//newline//&
         "The RUN_TYPE in the global section has to be set to NONE for FARMING."//newline//&
         "The different groups are executed in parallel. The jobs inside the same groups in series.",&
         repeats=.FALSE., required=.TRUE., supported_feature=.TRUE.,error=error)
    NULLIFY(keyword, print_key)

    CALL keyword_create(keyword, name="MASTER_SLAVE",&
         description="If a master-slave setup should be employed, in which one process is used to distribute the tasks. "//&
         "This is most useful to load-balance if not all jobs have the same length, "//&
         "and a lot of CPUs/groups are availabe.",&
         usage="MASTER_SLAVE",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,&
         supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="NGROUPS",variants=(/"NGROUP"/),&
         description="Gives the preferred number of working groups.",&
         usage="ngroups 4", type_of_var=integer_t, supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="GROUP_SIZE",&
         description="Gives the preferred size of a working group, "//&
         "groups will always be equal or larger than this size.",&
         usage="group_size 2", default_i_val=8, supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="GROUP_PARTITION",&
         description="gives the exact number of processors for each group.",&
         usage="group_partition  2 2 4 2 4 ", type_of_var=integer_t, n_var=-1, supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MAX_JOBS_PER_GROUP",&
         variants=(/"MAX_JOBS"/),&
         description="maximum number of jobs executed per group",&
         usage="max_step 4", default_i_val=HUGE(0)/32768, supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="CYCLE",&
         description="If farming should process all jobs in a cyclic way, stopping only if MAX_JOBS_PER_GROUP is exceeded." ,&
         usage="CYCLE",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,&
         supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="WAIT_TIME",&
         description="Time to wait [s] for a new task if no task is currently available, make this zero if no clock is available",&
         usage="WAIT_TIME 0.1",default_r_val=0.5_dp, supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    NULLIFY(sub_section)
    CALL section_create(sub_section,name="JOB",&
         description="description of the jobs to be executed",&
         repeats=.TRUE., required=.TRUE., supported_feature=.TRUE.,error=error)

    CALL keyword_create(keyword,name="DIRECTORY",&
         description="the directory in which the job should be executed",&
         usage="DIRECTORY /my/path",&
         required=.TRUE.,default_lc_val=".",supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(sub_section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword,name="INPUT_FILE_NAME",&
         description="the filename of the input file",&
         usage="INPUT_FILE_NAME my_input.inp",&
         required=.TRUE.,default_lc_val="input.inp",supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(sub_section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword,name="OUTPUT_FILE_NAME",&
         description="the filename of the output file, if not specified will use the project name in the &GLOBAL section.",&
         usage="OUTPUT_FILE_NAME my_input.inp",&
         required=.FALSE.,default_lc_val="",supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(sub_section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword,name="JOB_ID",&
         description="An ID used to indentify a job in DEPENDENCIES. "//&
         "JOB_IDs do not need to be unique, dependencies will be on all jobs with a given ID. "//&
         "If no JOB_ID is given, the index of the &JOB section in the input file will be used. ",&
         usage="JOB_ID 13", type_of_var=integer_t, supported_feature=.TRUE.,error=error)
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
         usage="DEPENDENCIES 13 1 7",type_of_var=integer_t, n_var=-1, supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(sub_section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL section_add_subsection(section, sub_section, error=error)
    CALL section_release(sub_section,error=error)

    CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
         description="controls the printing of FARMING specific output",&
         print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
         supported_feature=.TRUE.,error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL keyword_create(keyword, name="DO_RESTART",&
         description="Restart a farming job (and should pick up where the previous left off)",&
         usage="DO_RESTART",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,&
         supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="RESTART_FILE_NAME",&
         description="Name of the restart file to use for restarting a FARMING run. If not "//&
         "specified the name is determined from PROJECT name.",&
         usage="RESTART_FILE_NAME <FILENAME>", type_of_var=lchar_t,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL cp_print_key_section_create(print_key,"RESTART",&
         description="controls the printing of the restart for FARMING.",&
         print_level=low_print_level,add_last=add_last_numeric,filename="FARMING",&
         supported_feature=.TRUE.,error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

  END SUBROUTINE create_farming_section


! *****************************************************************************
!> \brief section for the chemical potential difference run
!> \param section will contain the cpot section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
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


! *****************************************************************************
!> \brief   creates the rs_pw_transfer section for use in the test section
!> \author  Joost VandeVondele
!> \date    2008-03-09
! *****************************************************************************
  SUBROUTINE create_rs_pw_transfer_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_rs_pw_transfer_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="RS_PW_TRANSFER",&
         description="Describes how to benchmark the rs_pw_transfer routines.",&
         n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
         error=error)

    NULLIFY(keyword)
    CALL keyword_create(keyword, name="GRID",&
         description="Specify the number of grid points (not all grid points are allowed)",&
         usage="GRID_DIMENSIONS 128 128 128", type_of_var=integer_t, n_var=3,&
         default_i_vals=(/128,128,128/), error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="HALO_SIZE",&
         description="number of grid points of the halo",&
         usage="HALO_SIZE 17",default_i_val=17,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="N_LOOP",&
         description="Number of rs_pw_transfers being timed",&
         usage="N_LOOP 100",default_i_val=10,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="RS2PW",&
         description="should the direction be rs2pw (pw2rs otherwise)",&
         usage="rs2pw TRUE",default_l_val=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    NULLIFY(subsection)
    CALL create_rsgrid_section(subsection,error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

  END SUBROUTINE create_rs_pw_transfer_section


! *****************************************************************************
!> \brief   creates the rs_pw_transfer section for use in the test section
!> \author  Joost VandeVondele
!> \date    2008-03-09
! *****************************************************************************
  SUBROUTINE create_pw_transfer_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_pw_transfer_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="PW_TRANSFER",&
         description="Benchmark and test the pw_transfer routines.",&
         n_keywords=1, n_subsections=0, repeats=.TRUE., required=.FALSE.,&
         error=error)

    NULLIFY(keyword)
    CALL keyword_create(keyword, name="GRID",&
         description="Specify the number of grid points (not all grid points are allowed)",&
         usage="GRID_DIMENSIONS 128 128 128", type_of_var=integer_t, n_var=3,&
         default_i_vals=(/128,128,128/), error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="N_LOOP",&
         description="Number of pw_transfers (backward&forward) being timed",&
         usage="N_LOOP 100",default_i_val=100,error=error)
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

    CALL keyword_create(keyword, name="PW_GRID_LAYOUT_ALL",&
         description="loop overal all PW_GRID_LAYOUTs that are compatible with a given number of CPUs ",&
         usage="PW_GRID_LAYOUT_ALL",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE., error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="DEBUG",&
         description="Do the FFT in debug mode in all cases",&
         usage="DEBUG",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE., error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="PW_GRID_LAYOUT",&
         description="Expert use only, leave the default..."//&
                     "Can be used to set the distribution for ray-distributed FFT.",&
         usage="PW_GRID_LAYOUT",&
         repeats=.FALSE.,required=.FALSE.,n_var=2,&
         default_i_vals=(/-1,-1/),&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="PW_GRID_BLOCKED",&
         description="Expert use only, leave the default..."//&
                     "Can be used to set the distribution in g-space for the pw grids and their FFT.",&
         usage="PW_GRID_BLOCKED FREE",&
         enum_c_vals=s2a("FREE","TRUE","FALSE"),&
         enum_desc=s2a("CP2K will select the optimal value","blocked","not blocked"),&
         enum_i_vals=(/do_pw_grid_blocked_free,do_pw_grid_blocked_true,do_pw_grid_blocked_false/),&
         default_i_val=do_pw_grid_blocked_false, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)


  END SUBROUTINE create_pw_transfer_section


! *****************************************************************************
!> \brief   creates the cp_fm_gemm section for use in the test section
!> \author  Joost VandeVondele
!> \date    2009-06-15
! *****************************************************************************
  SUBROUTINE create_cp_fm_gemm_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_cp_fm_gemm_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="CP_FM_GEMM",&
         description="Benchmark and test the cp_fm_gemm routines by multiplying C=A*B  ",&
         n_keywords=1, n_subsections=0, repeats=.TRUE., required=.FALSE.,&
         error=error)

    NULLIFY(keyword)
    CALL keyword_create(keyword, name="N_LOOP",&
         description="Number of cp_fm_gemm operations being timed (useful for small matrices).",&
         usage="N_LOOP 10",default_i_val=10,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="K",&
         description="Dimension 1 of C",&
         usage="A 1024",default_i_val=256,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="M",&
         description="Inner dimension M   ",&
         usage="A 1024",default_i_val=256,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="N",&
         description="Dimension 2 of C",&
         usage="A 1024",default_i_val=256,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="NROW_BLOCK",&
         description="block_size for rows",&
         usage="nrow_block 64",default_i_val=32,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="NCOL_BLOCK",&
         description="block_size for cols",&
         usage="nrow_block 64",default_i_val=32,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="ROW_MAJOR",&
         description="Use a row major blacs grid",&
         usage="ROW_MAJOR .FALSE.",default_l_val=.TRUE.,lone_keyword_l_val=.TRUE., error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="FORCE_BLOCKSIZE",&
         description="Forces the blocksize, even if this implies that a few processes might have no data",&
         usage="FORCE_BLOCKSIZE",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE., error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="GRID_2D",&
         description="Explicitly set the blacs 2D processor layout."//&
                     " If the product differs from the number of MPI ranks,"//&
                     " it is ignored and a default nearly square layout is used.", n_var=2,&
         usage="GRID_2D 64 16 ",default_i_vals=(/1, 1/),error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="TRANSA",&
         description="Transpose matrix A",&
         usage="TRANSA",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE., error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="TRANSB",&
         description="Transpose matrix B",&
         usage="TRANSB",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE., error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

  END SUBROUTINE create_cp_fm_gemm_section


! *****************************************************************************
!> \brief   creates the eigensolver section for use in the test section
!> \author  Joost VandeVondele
!> \date    2010-03-10
! *****************************************************************************
  SUBROUTINE create_eigensolver_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_eigensolver_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="EIGENSOLVER",&
         description="Benchmark and test the eigensolver routines.",&
         n_keywords=1, n_subsections=0, repeats=.TRUE., required=.FALSE.,&
         error=error)

    NULLIFY(keyword)
    CALL keyword_create(keyword, name="N",&
         description="Dimension of the square matrix",&
         usage="N 1024",default_i_val=256,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="N_LOOP",&
         description="Number of operations being timed (useful for small matrices).",&
         usage="N_LOOP 10",default_i_val=10,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="DIAG_METHOD",&
         description="Diagonalization strategy",&
         usage="DIAG_METHOD syevd",&
         enum_c_vals=s2a("syevd","syevx"),&
         enum_desc=s2a("(sca)lapacks syevd","(sca)lapacks syevx"),&
         enum_i_vals=(/ do_diag_syevd, do_diag_syevx/),&
         default_i_val=do_diag_syevd, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="EIGENVALUES",&
         description="number of eigenvalues to be computed (all=<0) ",&
         usage="EIGENVALUES 13",default_i_val=-1,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="INIT_METHOD",&
         description="Initialization approach",&
         usage="INIT_METHOD RANDOM",&
         enum_c_vals=s2a("random","read"),&
         enum_desc=s2a("use a random initial matrix", "read a matrix from file MATRIX"),&
         enum_i_vals=(/ do_mat_random, do_mat_read/),&
         default_i_val=do_mat_random, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

  END SUBROUTINE create_eigensolver_section


! *****************************************************************************
!> \brief   creates the cp_dbcsr section for use in the test section
!> \author  Urban Borstnik
!> \date    2010-02-08
! *****************************************************************************
  SUBROUTINE create_cp_dbcsr_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_cp_dbcsr_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="CP_DBCSR",&
         description="Benchmark and test the cp_dbcsr routines",&
         n_keywords=1, n_subsections=0, repeats=.TRUE., required=.FALSE.,&
         error=error)

    NULLIFY(keyword)
    CALL keyword_create(keyword, name="N_LOOP",&
         description="Number of operations being timed (useful for small matrices).",&
         usage="N_LOOP 10",default_i_val=10,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="DATA_TYPE",&
         description="Data type of the matrices",&
         usage="DATA_TYPE real_4",&
         default_i_val=dbcsr_type_real_8,&
         enum_c_vals=s2a("real_4", "real_8", "complex_4", "complex_8"),&
         enum_i_vals=(/dbcsr_type_real_4, dbcsr_type_real_8, &
                       dbcsr_type_complex_4, dbcsr_type_complex_8 /),&
         enum_desc=s2a(&
         "Real (Single Precision)",&
         "Real (Double Precision)",&
         "Complex (Single Precision)",&
         "Complex (Double Precision)"),&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="M",&
         description="Dimension 1 of C",&
         usage="A 1024",default_i_val=256,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="N",&
         description="Dimension 2 of C",&
         usage="A 1024",default_i_val=256,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="K",&
         description="Inner dimension M   ",&
         usage="A 1024",default_i_val=256,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="TRANSA",&
         description="Transpose matrix A",&
         usage="TRANSA",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE., error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="TRANSB",&
         description="Transpose matrix B",&
         usage="TRANSB",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE., error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="BS_M",&
         description="Row block sizes of C", n_var=-1,&
         usage="BS_M 1 13 2 5",default_i_vals=(/1, 13, 2, 15/),error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="BS_N",&
         description="Column block sizes of C", n_var=-1,&
         usage="BS_N 1 13 2 5",default_i_vals=(/1, 13, 2, 15/),error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="BS_K",&
         description="Block sizes of inner dimension", n_var=-1,&
         usage="BS_K 1 13 2 5",default_i_vals=(/1, 13, 2, 15/),error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="ATYPE",&
         description="Matrix A type",&
         usage="ATYPE N",default_c_val='N',error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="BTYPE",&
         description="Matrix B type",&
         usage="BTYPE N",default_c_val='N',error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="CTYPE",&
         description="Matrix C type",&
         usage="CTYPE N",default_c_val='N',error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="NPROC",&
         description="Number of processors to test", n_var=-1,&
         usage="NPROC 128 16 1",default_i_vals=(/0/),error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)


    CALL keyword_create(keyword, name="KEEPSPARSE",&
         description="Keep product sparse",&
         usage="KEEPSPARSE",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE., error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="ASPARSITY",&
         description="Sparsity of A matrix",&
         usage="ASPARSITY 70",default_r_val=0.0_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="BSPARSITY",&
         description="Sparsity of B matrix",&
         usage="ASPARSITY 80",default_r_val=0.0_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="CSPARSITY",&
         description="Sparsity of C matrix",&
         usage="ASPARSITY 90",default_r_val=0.0_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="ALPHA",&
         description="Multiplication factor",&
         usage="ALPHA 2.0",default_r_val=1.0_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="BETA",&
         description="Product premultiplication factor",&
         usage="BETA 1.0",default_r_val=0.0_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="FILTER_EPS",&
         description="Threshold for on-the-fly and final filtering.",&
         usage="FILTER_EPS 1.0",default_r_val=-1.0_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="DETERMINISTIC",&
         description="Attempt deterministic behavior (turn off for speed).",&
         usage="DETERMINISTIC",default_l_val=.TRUE.,lone_keyword_l_val=.TRUE., error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="ALWAYS_CHECKSUM",&
         description="perform a checksum after each multiplication",&
         usage="ALWAYS_CHECKSUM",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE., error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

  END SUBROUTINE create_cp_dbcsr_section


END MODULE input_cp2k
