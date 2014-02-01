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
  USE bibliography,                    ONLY: Frigo2005
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             empty_initial_variables,&
                                             parser_create,&
                                             parser_release
  USE cp_units,                        ONLY: cp_unit_set_create,&
                                             cp_unit_set_release,&
                                             cp_unit_set_type
  USE dbcsr_config,                    ONLY: &
       cudadrv_avoid_after_busy, cudadrv_binning_binsize, &
       cudadrv_binning_nbins, cudadrv_do_inhomogenous, &
       cudadrv_min_flop_process, cudadrv_min_flop_sort, &
       cudadrv_posterior_buffers, cudadrv_posterior_streams, &
       cudadrv_priority_buffers, cudadrv_priority_streams, &
       dbcsr_get_conf_combtypes, dbcsr_get_conf_max_ele_block, &
       dbcsr_get_conf_mm_driver, dbcsr_get_conf_mpi_mem, &
       dbcsr_get_conf_subcomm, dbcsr_get_conf_use_comm_thread, has_blas, &
       has_cuda, has_plasma, has_smm_gemm, mm_driver_blas, mm_driver_cuda, &
       mm_driver_matmul, mm_driver_plasma, mm_driver_smm, mm_name_blas, &
       mm_name_cuda, mm_name_matmul, mm_name_plasma, mm_name_smm, &
       multrec_limit
  USE dbcsr_error_handling
  USE dbcsr_types,                     ONLY: dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8
  USE f77_blas
  USE input_constants
  USE input_cp2k_atom,                 ONLY: create_atom_section
  USE input_cp2k_force_eval,           ONLY: create_force_eval_section
  USE input_cp2k_motion,               ONLY: create_motion_section
  USE input_cp2k_rsgrid,               ONLY: create_rsgrid_section
  USE input_cp2k_vib,                  ONLY: create_vib_section
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_optimize_basis,            ONLY: create_optimize_basis_section
  USE input_optimize_input,            ONLY: create_optimize_input_section
  USE input_parsing,                   ONLY: section_vals_parse
  USE input_section_types,             ONLY: &
       section_add_keyword, section_add_subsection, section_create, &
       section_release, section_type, section_vals_create, section_vals_get, &
       section_vals_get_subs_vals, section_vals_type
  USE input_val_types,                 ONLY: char_t,&
                                             integer_t,&
                                             lchar_t,&
                                             logical_t,&
                                             real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: newline,&
                                             s2a
  USE swarm_input,                     ONLY: create_swarm_section
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k'

  PUBLIC :: create_cp2k_root_section, parsed_cp2k_input,&
       create_global_section, create_cp2k_input_reading, empty_initial_variables

CONTAINS

! *****************************************************************************
!> \brief reads the cp2k input from the given filepath and returns a section_vals
!>      containing the input
!> \param file_path path where the input should be read
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  FUNCTION create_cp2k_input_reading(file_path,initial_variables,para_env,error) RESULT(res)
    CHARACTER(len=*), INTENT(in)             :: file_path
    CHARACTER(len=*), DIMENSION(:, :)        :: initial_variables
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

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(res)
    NULLIFY(cpparser, input_structure, default_units)
    CALL create_cp2k_root_section(input_structure,error=error)
    CALL section_vals_create(res,input_structure, error=error)
    CALL section_release(input_structure,error=error)
    CALL parser_create(cpparser,initial_variables=initial_variables,file_name=file_path, &
         para_env=para_env, error=error)
    CALL cp_unit_set_create(default_units, "OUTPUT",error=error)
    CALL section_vals_parse(res,cpparser,root_section=.FALSE.,&
         default_units=default_units,error=error)
    CALL cp_unit_set_release(default_units,error=error)
    CALL parser_release(cpparser,error=error)
    CALL timestop(handle)
  END FUNCTION create_cp2k_input_reading

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
!> \brief utility function to ease the transition to the new input.
!>      returns true if the new input was parsed
!> \param input_file the parsed input file
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
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

! *****************************************************************************
!> \brief section to hold global settings for the whole program
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE create_global_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_global_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, sub_section

    failure=.FALSE.
    NULLIFY(print_key)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="GLOBAL",&
         description="Section with general information regarding which kind "//&
         "of simulation to perform an parameters for the whole PROGRAM",&
         n_keywords=7, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
         supported_feature=.TRUE.,error=error)

    NULLIFY(keyword)
    CALL keyword_create(keyword, name="BLACS_GRID",&
         description="how to distribute the processors on the 2d grid needed "//&
         "by BLACS (and thus SCALAPACK)", usage="BLACS_GRID SQUARE",&
         default_i_val=BLACS_GRID_SQUARE,enum_c_vals=s2a("SQUARE","ROW","COLUMN"),&
         enum_desc=s2a("Distribution by matrix blocks", "Distribution by matrix rows",&
         "Distribution by matrix columns"), &
         enum_i_vals=(/BLACS_GRID_SQUARE,BLACS_GRID_ROW,BLACS_GRID_COL/), error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="BLACS_REPEATABLE",&
         description="Use a topology for BLACS collectives that is guaranteed to be repeatable "//&
                     "on homegeneous architectures",&
         usage="BLACS_REPEATABLE",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="PREFERRED_DIAG_LIBRARY",&
         description="Specifies the DIAGONALIZATION library. If not availabe ....",&
         usage="PREFERRED_DIAG_LIBRARY ELPA",&
         default_i_val=do_diag_sl, &
         enum_i_vals=(/do_diag_sl,do_diag_sl2,do_diag_elpa/),&
         enum_c_vals=s2a("SL","SL2","ELPA"),&
         enum_desc=s2a("Standard scalapack: syevd","Scalapack 2.0: syevr","ELPA"),&
         supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)


    CALL keyword_create(keyword, name="PREFERRED_FFT_LIBRARY",&
         description="Specifies the FFT library which should be preferred. "//&
                     "If it is not available, use FFTW3 if this is linked in, if FFTW3 is not available use FFTSG. "//&
                     "Improved performance with FFTW3 can be obtained specifying a proper value for FFTW_PLAN_TYPE. "//&
                     "Contrary to earlier CP2K versions, all libraries will result in the same grids, "//&
                     "i.e. the subset of grids which all FFT libraries can transform. "//&
                     "See EXTENDED_FFT_LENGTHS if larger FFTs or grids that more precisely match a given cutoff are needed, "//&
                     "or older results need to be reproduced. "//&
                     "FFTW3 is often (close to) optimal, and well tested with CP2K.",&
         usage="PREFERRED_FFT_LIBRARY FFTW3",&
         citations=(/Frigo2005/),&
         default_i_val=do_fft_fftw3, &
         enum_i_vals=(/do_fft_sg,do_fft_fftw3,do_fft_fftw3/),&
         enum_c_vals=s2a("FFTSG","FFTW3","FFTW"),&
         enum_desc=s2a("Stefan Goedecker''s FFT (FFTSG), always available,"//&
                       "will be used in case a FFT lib is specified and not available",&
                       "a fast portable FFT library. Recommended."//&
                       "See also the FFTW_PLAN_TYPE, and FFTW_WISDOM_FILE_NAME keywords.",&
                       "Same as FFTW3 (for compatability with CP2K 2.3)"),&
         supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="FFTW_WISDOM_FILE_NAME",&
         description="The name of the file that contains wisdom (pre-planned FFTs) for use with FFTW3. "//&
                     "Using wisdom can significantly speed up the FFTs (see the FFTW homepage for details). "//&
                     "Note that wisdom is not transferable between different computer (architectures). "//&
                     "Wisdom can be generated using the fftw-wisdom tool that is part of the fftw installation. "//&
                     "cp2k/tools/cp2k-wisdom is a script that contains some additional info, and can help "//&
                     "to generate a useful default for /etc/fftw/wisdom or particular values for a given simulation.",&
         usage="FFTW_WISDOM_FILE_NAME wisdom.dat", default_lc_val="/etc/fftw/wisdom", supported_feature=.TRUE.,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

   CALL keyword_create(keyword, name="FFTW_PLAN_TYPE",&
         description="FFTW can have improved performance if it is allowed to plan with "//&
                     "explicit measurements which strategy is best for a given FFT. "//&
                     "While a plan based on measurements is generally faster, "//&
                     "differences in machine load will lead to different plans for the same input file, "//&
                     "and thus numerics for the FFTs will be slightly different from run to run."//&
                     "PATIENT planning is recommended for long ab initio MD runs.",&
         usage="FFTW_PLAN_TYPE PATIENT",&
         citations=(/Frigo2005/),&
         default_i_val=fftw_plan_estimate, &
         enum_i_vals=(/fftw_plan_estimate,fftw_plan_measure,fftw_plan_patient,fftw_plan_exhaustive/), &
         enum_c_vals=s2a("ESTIMATE",&
                         "MEASURE",&
                         "PATIENT",&
                         "EXHAUSTIVE"),&
         enum_desc=s2a("Quick estimate, no runtime measurements.",&
                       "Quick measurement, somewhat faster FFTs.",&
                       "Measurements trying a wider range of possibilities.",&
                       "Measurements trying all possibilities - use with caution."),&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="EXTENDED_FFT_LENGTHS",&
         description="Use fft library specific values for the allows number of points in FFTs. "//&
                     "The default is to use the internal FFT lengths. For external fft libraries this may "//&
                     "create an error at the external library level, because the length provided by cp2k is "//&
                     "not supported by the external library. In this case switch on this keyword "//&
                     "to obtain, with certain fft libraries, lengths matching the external fft library lengths, or "//&
                     "larger allowed grids, or grids that more precisely match a given cutoff. "//&
                     "IMPORTANT NOTE: in this case, the actual grids used in CP2K depends on the FFT library. "//&
                     "A change of FFT library must therefore be considered equivalent to a change of basis, "//&
                     "which implies a change of total energy. ",&
         usage="EXTENDED_FFT_LENGTHS",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword,name="FFT_POOL_SCRATCH_LIMIT",&
         description="Limits the memory usage of the FFT scratch pool, potentially reducing efficiency a bit",&
         usage="FFT_POOL_SCRATCH_LIMIT <INTEGER>",default_i_val=15,error=error)
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
         high_print_level,debug_print_level/),supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
                   
    CALL keyword_create(keyword, name="PROGRAM_NAME",&
         variants=(/"PROGRAM"/),&
         description="Which program should be run",&
         usage="PROGRAM_NAME <STRING>",&
         enum_c_vals=s2a("ATOM","FARMING","TEST","CP2K","OPTIMIZE_INPUT","OPTIMIZE_BASIS", "TMC","MC_ANALYSIS","SWARM"),&
         enum_desc=s2a("Runs single atom calculations",&
                       "Runs N independent jobs in a single run",&
                       "Do some benchmarking and testing",&
                       "Runs one of the CP2K package",&
                       "A tool do optimize parameters in a CP2K input", &
                       "A tool to create a MOLOPT or ADMM basis for a given set"//&
                       " of training structures",&
                       "Runs Tree Monte Carlo algorithm using additional input file(s)", &
                       "Runs (Tree) Monte Carlo trajectory file analysis",&
                       "Runs swarm based calculation"),&
         enum_i_vals=(/do_atom, do_farming, do_test, do_cp2k, do_optimize_input,&
                       do_opt_basis, do_tree_mc, do_tree_mc_ana, do_swarm/),&
         default_i_val=do_cp2k, supported_feature=.TRUE., error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="PROJECT_NAME",&
         variants=(/"PROJECT"/),&
         description="Name of the project (used to build the name of the "//&
         "trajectory, and other files generated by the program)",&
         usage="PROJECT_NAME <STRING>",&
         default_c_val="PROJECT",supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="OUTPUT_FILE_NAME",&
         description="Name of the output file. "//&
         "Relevant only if automatically started (through farming for example). "//&
         "If empty uses the project name as basis for it.",&
         usage="OUTPUT_FILE_NAME {filename}",default_lc_val="", supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="RUN_TYPE",&
         description="Type of run that you want to perform Geometry "//&
         "optimization, md, montecarlo,...",&
         usage="RUN_TYPE MD",&
         default_i_val=energy_force_run,&
         enum_c_vals=s2a("NONE", "ENERGY", "ENERGY_FORCE", "MD", "GEO_OPT",&
         "MC", "SPECTRA", "DEBUG", "BSSE", "LR", "PINT", "VIBRATIONAL_ANALYSIS",&
         "BAND", "CELL_OPT", "WFN_OPT", "WAVEFUNCTION_OPTIMIZATION",&
         "MOLECULAR_DYNAMICS", "GEOMETRY_OPTIMIZATION", "MONTECARLO",&
         "ELECTRONIC_SPECTRA", "LINEAR_RESPONSE", "NORMAL_MODES","RT_PROPAGATION",&
         "EHRENFEST_DYN","TAMC","TMC","DRIVER"),&
         enum_i_vals=(/ none_run, energy_run, energy_force_run, mol_dyn_run,&
          geo_opt_run, mon_car_run, electronic_spectra_run, debug_run,&
          bsse_run, linear_response_run, pint_run, vib_anal,do_band,&
          cell_opt_run, energy_run, energy_run, mol_dyn_run, geo_opt_run,&
          mon_car_run, electronic_spectra_run, linear_response_run,&
          vib_anal,real_time_propagation,ehrenfest,do_tamc, tree_mc_run,driver_run/),&
         enum_desc=s2a("Perform no tasks", "Computes energy","Computes energy and forces",&
         "Molecular Dynamics","Geometry Optimization","Monte Carlo", "Computes absorption Spectra",&
         "Performs a Debug analysis","Basis set superposition error","Linear Response",&
         "Path integral","Vibrational analysis","Band methods","Cell optimization",&
         "Alias for ENERGY","Alias for ENERGY","Alias for MD","Alias for GEO_OPT",&
         "Alias for MC","Alias for SPECTRA","Alias for LR","Alias for VIBRATIONAL_ANALYSIS",&
         "Real Time propagation run (fixed ionic positions)",&
         "Ehrenfest dynamics (using real time propagation of the wavefunction)",&
         "TAMC", "TreeMonteCarlo (TMC), a pre sampling MC algorithma","i-PI driver mode"),&
         supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="WALLTIME",&
          variants=(/"WALLTI"/),&
          description="Maximum execution time for this run. Time in seconds or in HH:MM:SS.",&
          usage="WALLTIME {real} or {HH:MM:SS}",  default_lc_val="", &
          supported_feature=.TRUE.,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="ECHO_INPUT",&
         description="If the input should be echoed to the output with all the "//&
         "defaults made explicit",&
         usage="ECHO_INPUT NO",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,&
         supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="ECHO_ALL_HOSTS",&
         description="Echo a list of hostname and pid for all MPI processes.",&
         usage="ECHO_ALL_HOSTS NO",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="TRACE",&
         description="If a debug trace of the execution of the program should be written ",&
         usage="TRACE",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    
    CALL keyword_create(keyword, name="TRACE_MASTER",&
         description="For parallel TRACEd runs: only the master node writes output.",&
         usage="TRACE_MASTER",&
         default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="TRACE_MAX",&
         description="Limit the total number a given subroutine is printed in the trace. Accounting is not influenced.",&
         usage="TRACE_MAX 100",default_i_val=HUGE(0),error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    
    CALL keyword_create(keyword, name="TRACE_ROUTINES",&
         description="A list of routines to trace. If left empty all routines are traced. Accounting is not influenced.",&
         usage="TRACE_ROUTINES {routine_name1} {routine_name2} ...", type_of_var=char_t,&
         n_var=-1, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    
    CALL keyword_create(keyword, name="FLUSH_SHOULD_FLUSH",&
         description="Flush output regularly, enabling this option might degrade performance significantly on certain machines.",&
         usage="FLUSH_SHOULD_FLUSH",&
         default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    
    CALL keyword_create(keyword, name="CALLGRAPH",&
    description="At the end of the run write a callgraph to file, "//&
         "which contains detailed timing informations. "//&
         "This callgraph can be viewed e.g. with the open-source program kcachegrind.",&
         usage="CALLGRAPH <NONE|MASTER|ALL>",&
         default_i_val=CALLGRAPH_NONE, lone_keyword_i_val=CALLGRAPH_MASTER,&
         enum_c_vals=s2a("NONE","MASTER","ALL"),&
         enum_desc=s2a("No callgraph gets written",&
         "Only the master process writes his callgraph",&
         "All processes write their callgraph (into a separate files)."), &
         enum_i_vals=(/CALLGRAPH_NONE, CALLGRAPH_MASTER, CALLGRAPH_ALL/), error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    
    CALL keyword_create(keyword, name="CALLGRAPH_FILE_NAME",&
         description="Name of the callgraph file, which is writte a the end of the run. "//&
         "If not specified the project name will be used as filename.",&
         usage="CALLGRAPH_FILE_NAME {filename}",default_lc_val="", supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    
    CALL keyword_create(keyword,name="SEED",&
         description="Initial seed for the global (pseudo)random number "//&
         "generator to create a stream of normally Gaussian "//&
         "distributed random numbers.",&
         usage="SEED <INTEGER>",default_i_val=2000,supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="SAVE_MEM",&
         description="Some sections of the input structure are deallocated when not needed,"//&
         " and reallocated only when used. This reduces the required maximum memory  ",&
         usage="SAVE_MEM",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="ENABLE_UNSUPPORTED_FEATURES",&
         description="This keywords enables the usage of unsupported features "//&
         "in a release version. It affects ONLY release versions of CP2K (no effects on"//&
         " development versions).",&
         usage="ENABLE_UNSUPPORTED_FEATURES <LOGICAL>",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL cp_print_key_section_create(print_key, "TIMINGS", description=&
         "Controls the printing of the timing report at the end of CP2K execution", &
         print_level=silent_print_level,filename="__STD_OUT__",&
         supported_feature=.TRUE.,error=error)
    CALL keyword_create(keyword,name="THRESHOLD",&
         description="Specify % of CPUTIME above which the contribution will be inserted in the"//&
         " final timing report (e.g. 0.02 = 2%)",&
         usage="THRESHOLD <real>",&
         default_r_val=0.02_dp, supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword,name="SORT_BY_SELF_TIME",&
         description="Sort the final timing report by the average self (exclusive) time instead of the "//&
         "total (inclusive) time of a routine",&
         usage="SORT_BY_SELF_TIME on",& 
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key, "REFERENCES", description=&
         "Controls the printing of the references relevant to the calculations performed", &
         print_level=silent_print_level,filename="__STD_OUT__",&
         supported_feature=.TRUE.,error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
         description="controls the printing of  initialization controlled by the global section",&
            print_level=silent_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            supported_feature=.TRUE.,error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key, "PRINT", description=&
         "controls the printing of physical and mathematical constants", &
         print_level=medium_print_level,filename="__STD_OUT__",&
         supported_feature=.TRUE.,error=error)

    CALL keyword_create(keyword,"BASIC_DATA_TYPES",&
         description="Controls the printing of the basic data types.",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword,"physcon",&
         description="if the printkey is active prints the physical constants",&
         default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword,"spherical_harmonics",&
         description="if the printkey is active prints the spherical harmonics",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword,"RNG_MATRICES",&
                        description="Prints the transformation matrices used by the "//&
                                    " random number generator",&
                        default_l_val=.FALSE.,&
                        lone_keyword_l_val=.TRUE.,&
                        supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword,"RNG_CHECK",&
                        description="Performs a check of the global (pseudo)random "//&
                                    "number generator (RNG) and prints the result",&
                        default_l_val=.FALSE.,&
                        lone_keyword_l_val=.TRUE.,&
                        supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword,"GLOBAL_GAUSSIAN_RNG",&
                        description="Prints the initial status of the global Gaussian "//&
                                    "(pseudo)random number stream which is mostly used for "//&
                                    "the velocity initialization",&
                        default_l_val=.FALSE.,&
                        lone_keyword_l_val=.TRUE.,&
                        supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword,"HPM",&
         description="if the keyword is active controls the printing of the hpm info (if linked as"//&
         " external library, on IBM only)",&
         default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)
! input for CUDA support
    NULLIFY(sub_section)
    CALL section_create(sub_section,name="CUDA",&
         description="For CUDA support related input",&
         repeats=.FALSE., required=.FALSE., supported_feature=.TRUE.,error=error)
    CALL keyword_create(keyword, name="MEMORY",&
         description="Allocates this much device memory (in kB)",&
         usage="MEMORY 10",default_i_val=-1,error=error)
    CALL section_add_keyword(sub_section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL section_add_subsection(section,sub_section,error=error)
    CALL section_release(sub_section,error=error)
    ! Machine arch
    CALL create_march_section(sub_section,error)
    CALL section_add_subsection(section,sub_section,error=error)
    CALL section_release(sub_section,error=error)
    ! DBCSR options
    CALL create_dbcsr_section(sub_section,error)
    CALL section_add_subsection(section,sub_section,error=error)
    CALL section_release(sub_section,error=error)
 END SUBROUTINE create_global_section

! *****************************************************************************
!> \brief   Creates the machine architecture section
!> \author  Christiane Pousa Ribeiro
!> \date    2011-08-15
! *****************************************************************************
  SUBROUTINE create_march_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_march_section', &
      routineP = moduleN//':'//routineN
    CHARACTER(len=32), PARAMETER :: has = "", &
      hasnot = " (not available in this binary)"

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="MACHINE_ARCH",&
         description="Configuration options for the machine architecture.",&
         n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
         error=error)
    !
    NULLIFY(keyword)
    CALL keyword_create(keyword, name="PRINT_RESUME",&
         description="Print the compute node architecture number of components.", &
         usage="PRINT_RESUME TRUE",&
         default_l_val=.FALSE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    !
    NULLIFY(keyword)
    CALL keyword_create(keyword, name="PRINT_FULL",&
         description="Print full machine architecture", &
         usage="PRINT_FULL TRUE",&
         default_l_val=.FALSE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

     !
    NULLIFY(keyword)
    CALL keyword_create(keyword, name="PRINT_BRANCH",&
         description="Print machine architecture component organization.", &
         usage="PRINT_BRANCH TRUE",&
         default_l_val=.FALSE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
      !
    NULLIFY(keyword)
    CALL keyword_create(keyword, name="PRINT_THREAD",&
         description="Print the threads binding in the machine.", &
         usage="PRINT_THREAD TRUE",&
         default_l_val=.FALSE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
      !
    NULLIFY(keyword)
    CALL keyword_create(keyword, name="PRINT_THREAD_CUR",&
         description="Print the current threads binding in the machine.", &
         usage="PRINT_THREAD_CUR TRUE",&
         default_l_val=.FALSE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

   !
    NULLIFY(keyword)
    CALL keyword_create(keyword, name="PRINT_PROC",&
         description="Print the process binding in the machine.", &
         usage="PRINT_PROC TRUE",&
         default_l_val=.FALSE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
   !
    NULLIFY(keyword)
    CALL keyword_create(keyword, name="SCHED_THREAD",&
         description="Enable thread scheduling on the compute node.", &
         usage="SCHED_THREAD type (N=none, L=linear, I=interleaved, "//&
          "G=group,M=manual,D=default)",&
         default_c_val='D',error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
  !
    NULLIFY(keyword)
    CALL keyword_create(keyword, name="SCHED_MPI",&
         description="Enable process scheduling on the compute node.", &
         usage="SCHED_MPI type (N=none, L=linear, I=interleaved, G=group,"//&
         " M=manual, D=Default)",&
         default_c_val='D',error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
   !
    NULLIFY(keyword)
    CALL keyword_create(keyword, name="MEM_POL",&
         description="Enable memory binding on the compute node.", &
         usage="MEM_POL type (D=default, L=local, I=interleave, M=manual )",&
         default_c_val='X',error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
   !
   NULLIFY(keyword)
   CALL keyword_create(keyword, name="PROCESS",&
        description="Core ID for each MPI process of a compute node. It will be"//&
                    " used to place the MPI process.", n_var=-1,&
        usage="PROCESS 0 1 2 NMPI",default_i_vals=(/0/),error=error)
   CALL section_add_keyword(section,keyword,error=error)
   CALL keyword_release(keyword,error=error)
   !
   NULLIFY(keyword)
   CALL keyword_create(keyword, name="THREAD",&
        description="Core ID for each thread of a compute node. It will be"//&
                     " used to place the thread.", n_var=-1,&
        usage="THREAD 0 1 2 NThreads",default_i_vals=(/0/),error=error)
   CALL section_add_keyword(section,keyword,error=error)
   CALL keyword_release(keyword,error=error)
   !
   NULLIFY(keyword)
   CALL keyword_create(keyword, name="MEMORY",&
        description="The memory banks ID for each MPI process. It will be"//&
                    "use to set memory bank of the MPI process", n_var=-1,&
        usage="MEMORY 0 1 2 NMPI",default_i_vals=(/0/),error=error)
   CALL section_add_keyword(section,keyword,error=error)
   CALL keyword_release(keyword,error=error)
   !
   NULLIFY(keyword)
   CALL keyword_create(keyword, name="MPI_REORDERING",&
        description="Apply a global MPI reordering for the run", n_var=-1,&
        usage="MPI_REORDERING type (D=none H=hilbert P=peano N=snake curve " //&
              "R=round-robin, F=hilbert-peano, S=switch pairs, C=cannon, O=own ordering) ",&
        default_c_val='N',error=error)
   CALL section_add_keyword(section,keyword,error=error)
   CALL keyword_release(keyword,error=error)
   

  END SUBROUTINE create_march_section

! *****************************************************************************
!> \brief   Creates the dbcsr section for configuring DBCSR
!> \author  Urban Borstnik
!> \date    2011-04-05
! *****************************************************************************
  SUBROUTINE create_dbcsr_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_dbcsr_section', &
      routineP = moduleN//':'//routineN
    CHARACTER(len=32), PARAMETER             :: has = "", &
                                                hasnot = " (not available)"

    CHARACTER(len=64)                        :: a_blas, a_cuda, a_plasma, &
                                                a_smm
    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="DBCSR",&
         description="Configuration options for the DBCSR library.",&
         n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
         error=error)
    !
    IF (has_blas) THEN
       a_blas = has
    ELSE
       a_blas = hasnot
    ENDIF
    IF (has_smm_gemm) THEN
       a_smm = has
    ELSE
       a_smm = " (only BLAS wrapper available in this binary)"
    ENDIF
    IF (has_plasma) THEN
       a_plasma = has
    ELSE
       a_plasma = hasnot
    ENDIF
    IF (has_cuda) THEN
       a_cuda = has
    ELSE
       a_cuda = hasnot
    ENDIF
    !
    NULLIFY(keyword)
    CALL keyword_create(keyword, name="subcommunicators",&
         description="Use MPI subcommunicators for transfers "// &
         "that are limited to process grid rows and columns.",&
         usage="subcommunicators T",&
         default_l_val=dbcsr_get_conf_subcomm(),error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    !
    CALL keyword_create(keyword, name="combined_types",&
         description="Combine matrix index and data into "// &
         "a single MPI derived data type for communication (avoid using).",&
         usage="combined_types T",&
         default_l_val=dbcsr_get_conf_combtypes(),error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    !
    CALL keyword_create(keyword, name="use_MPI_allocation",&
         description="Allow use of MPI-allocated memory "//&
         "for potentially faster network communication.",&
         usage="use_MPI_allocation T",&
         default_l_val=dbcsr_get_conf_mpi_mem(),error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    !
    CALL keyword_create(keyword, name="mm_stack_size",&
         description="Size of multiplication parameter stack."&
         //" A negative value leaves the decision up to DBCSR.",&
         usage="mm_stack_size 1000",&
         default_i_val=-1,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    !
    CALL keyword_create(keyword, name="mm_driver",&
         description="Select which routines to use "//&
         "for matrix block multiplications.",&
         usage="mm_driver blas",&
         default_i_val=dbcsr_get_conf_mm_driver(),&
         enum_c_vals=s2a(mm_name_blas,&
         mm_name_matmul,&
         mm_name_smm,&
         mm_name_plasma,&
         mm_name_cuda),&
         enum_i_vals=(/&
         mm_driver_blas,&
         mm_driver_matmul,&
         mm_driver_smm,&
         mm_driver_plasma,&
         mm_driver_cuda/),&
         enum_desc=s2a(&
         TRIM("BLAS"//a_blas),&
         TRIM("Fortran MATMUL"),&
         TRIM("libsmm: a library for small matrix multiplies"//a_smm),&
         TRIM("PLASMA"//a_plasma),&
         TRIM("CUDA"//a_cuda)&
         ),&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    !
    CALL keyword_create(keyword, name="n_size_mnk_stacks",&
         description="Number of stacks to use for distinct atomic sizes" &
         // " (e.g., 2 for a system of mostly waters). "&
         //"A negative value leaves the decision up to DBCSR.",&
         usage="n_size_mnk_stacks 2",&
         default_i_val=-1,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    !
    CALL keyword_create(keyword, name="use_comm_thread",&
         description="During multiplication, use a thread to periodically poll" &
         // " MPI to progress outstanding message completions.  This is" &
         // " beneficial on systems without a DMA-capable network adapter" &
         // " e.g. Cray XE6.",&
         usage="use_comm_thread T",&
         default_l_val=dbcsr_get_conf_use_comm_thread(),error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    !
    CALL keyword_create(keyword, name="MAX_ELEMENTS_PER_BLOCK",&
         description="Default block size for turning dense matrices in blocked ones",&
         usage="MAX_ELEMENTS_PER_BLOCK 32",&
         default_i_val=dbcsr_get_conf_max_ele_block(),error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="comm_thread_load",&
         description="If a communications thread is used, specify how much " &
         // "multiplication workload (%) the thread should perform in " &
         // "addition to communication tasks. "&
         // "A negative value leaves the decision up to DBCSR.",&
         usage="comm_thread_load 50",&
         default_i_val=-1, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="multrec_limit",&
         description="Recursion limit of cache oblivious multrec algorithm.",&
         default_i_val=multrec_limit,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    !---------------------------------------------------------------------------
    NULLIFY(subsection)
    CALL section_create(subsection,name="CUDA",&
         description="Configuration options for the CUDA-Driver.",&
         n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
         error=error)

    CALL keyword_create(keyword, name="priority_buffers",&
         description="Number of transfer-buffers associated with high priority streams.",&
         default_i_val=cudadrv_priority_buffers,error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="posterior_buffers",&
         description="Number of transfer-buffers associated with low priority streams.",&
         default_i_val=cudadrv_posterior_buffers,error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="priority_streams",&
         description="Number of cuda streams created with high priority.",&
         default_i_val=cudadrv_priority_streams,error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="posterior_streams",&
         description="Number of cuda streams created with low priority.",&
         default_i_val=cudadrv_posterior_streams,error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="avoid_after_busy",&
         description="If enabled, stacks are not processed by the cuda-driver "&
         //"after it has signaled congestion during a round of flushing. "&
         //"For the next round of flusing the driver is used again.",&
         default_l_val=cudadrv_avoid_after_busy,error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="min_flop_process",&
         description="Only process stacks with more than the given number of "&
         //"floating-point operations per stack-entry (2*m*n*k).",&
         default_i_val=cudadrv_min_flop_process,error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="min_flop_sort",&
         description="Only sort stacks with more than the given number of "&
         //"floating-point operations per stack-entry (2*m*n*k). "&
         //"Alternatively, the stacks are roughly ordered through a "&
         //"binning-scheme by Peter Messmer.",&
         default_i_val=cudadrv_min_flop_sort,error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="process_inhomogenous",&
         description="If enabled, inhomogenous stacks are also processed by the cuda driver.",&
         default_l_val=cudadrv_do_inhomogenous,error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="binning_nbins",&
         description="Number of bins used when ordering "&
         //"the stacks with the binning-scheme.",&
         default_i_val=cudadrv_binning_nbins,error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="binning_binsize",&
         description="Size of bins used when ordering "&
         //"the stacks with the binning-scheme.",&
         default_i_val=cudadrv_binning_binsize,error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL section_add_subsection(section, subsection,error=error)
    CALL section_release(subsection, error=error)
    !
  END SUBROUTINE create_dbcsr_section


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

END MODULE input_cp2k
