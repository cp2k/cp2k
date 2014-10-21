!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief builds the global input section for cp2k
!> \par History
!>      06.2004 created [fawzi]
!>      03.2014 moved to seperate file [Ole Schuett]
!> \author fawzi
! *****************************************************************************
MODULE input_cp2k_global
  USE bibliography,                    ONLY: Borstnik2014,&
                                             Frigo2005
  USE cp_blacs_env,                    ONLY: BLACS_GRID_COL,&
                                             BLACS_GRID_ROW,&
                                             BLACS_GRID_SQUARE
  USE cp_dbcsr_interface,              ONLY: &
       accdrv_avoid_after_busy, accdrv_binning_binsize, accdrv_binning_nbins, &
       accdrv_do_inhomogenous, accdrv_min_flop_process, accdrv_min_flop_sort, &
       accdrv_posterior_buffers, accdrv_posterior_streams, &
       accdrv_priority_buffers, accdrv_priority_streams, &
       dbcsr_get_conf_combtypes, dbcsr_get_conf_max_ele_block, &
       dbcsr_get_conf_mm_driver, dbcsr_get_conf_mpi_mem, &
       dbcsr_get_conf_subcomm, dbcsr_get_conf_use_comm_thread, mm_driver_acc, &
       mm_driver_blas, mm_driver_matmul, mm_driver_smm, mm_name_acc, &
       mm_name_blas, mm_name_matmul, mm_name_smm, multrec_limit
  USE cp_output_handling,              ONLY: add_last_numeric,&
                                             cp_print_key_section_create,&
                                             debug_print_level,&
                                             high_print_level,&
                                             low_print_level,&
                                             medium_print_level,&
                                             silent_print_level
  USE f77_blas
  USE input_constants,                 ONLY: &
       bsse_run, callgraph_all, callgraph_master, callgraph_none, &
       cell_opt_run, debug_run, do_atom, do_band, do_cp2k, do_dbcsr, &
       do_diag_elpa, do_diag_sl, do_diag_sl2, do_farming, do_fft_fftw3, &
       do_fft_sg, do_opt_basis, do_optimize_input, do_pdgemm, do_swarm, &
       do_tamc, do_test, do_tree_mc, do_tree_mc_ana, driver_run, ehrenfest, &
       electronic_spectra_run, energy_force_run, energy_run, &
       fftw_plan_estimate, fftw_plan_exhaustive, fftw_plan_measure, &
       fftw_plan_patient, gaussian, geo_opt_run, linear_response_run, &
       mol_dyn_run, mon_car_run, none_run, pint_run, real_time_propagation, &
       tree_mc_run, vib_anal
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: char_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_global'

  PUBLIC :: create_global_section

CONTAINS

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
         "ELECTRONIC_SPECTRA", "LINEAR_RESPONSE", "NORMAL_MODES", "RT_PROPAGATION",&
         "EHRENFEST_DYN", "TAMC", "TMC", "DRIVER"),&
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
         "Temperature Accelerated Monte Carlo (TAMC)",&
         "Tree Monte Carlo (TMC), a pre-sampling MC algorithm",&
         "i-PI driver mode"),&
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
    NULLIFY(sub_section)
    ! FM section
    CALL create_fm_section(sub_section,error)
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
!> \param section ...
!> \param error ...
!> \date    2011-08-15
!> \author  Christiane Pousa Ribeiro
! *****************************************************************************
  SUBROUTINE create_march_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_march_section', &
      routineP = moduleN//':'//routineN

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
!> \param section ...
!> \param error ...
!> \date    2011-04-05
!> \author  Urban Borstnik
! *****************************************************************************
  SUBROUTINE create_dbcsr_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_dbcsr_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="DBCSR",&
         description="Configuration options for the DBCSR library.",&
         n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
         citations=(/Borstnik2014/), error=error)

    NULLIFY (keyword)

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
         enum_c_vals=s2a(mm_name_blas,mm_name_matmul,mm_name_smm,&
                         mm_name_acc),&
         enum_i_vals=(/mm_driver_blas,mm_driver_matmul,mm_driver_smm,&
                       mm_driver_acc/),&
         enum_desc=s2a("BLAS (requires the BLAS library at link time)",&
                       "Fortran MATMUL",&
                       "Library optimised for Small Matrix Multiplies "//&
                       "(requires the SMM library at link time)",&
                       "ACC (requires an accelerator backend)"),&
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
    CALL section_create(subsection,name="ACC",&
         description="Configuration options for the ACC-Driver.",&
         n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
         error=error)

    CALL keyword_create(keyword, name="priority_buffers",&
         description="Number of transfer-buffers associated with high priority streams.",&
         default_i_val=accdrv_priority_buffers,error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="posterior_buffers",&
         description="Number of transfer-buffers associated with low priority streams.",&
         default_i_val=accdrv_posterior_buffers,error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="priority_streams",&
         description="Number of acc streams created with high priority.",&
         default_i_val=accdrv_priority_streams,error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="posterior_streams",&
         description="Number of acc streams created with low priority.",&
         default_i_val=accdrv_posterior_streams,error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="avoid_after_busy",&
         description="If enabled, stacks are not processed by the acc-driver "&
         //"after it has signaled congestion during a round of flushing. "&
         //"For the next round of flusing the driver is used again.",&
         default_l_val=accdrv_avoid_after_busy,error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="min_flop_process",&
         description="Only process stacks with more than the given number of "&
         //"floating-point operations per stack-entry (2*m*n*k).",&
         default_i_val=accdrv_min_flop_process,error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="min_flop_sort",&
         description="Only sort stacks with more than the given number of "&
         //"floating-point operations per stack-entry (2*m*n*k). "&
         //"Alternatively, the stacks are roughly ordered through a "&
         //"binning-scheme by Peter Messmer.",&
         default_i_val=accdrv_min_flop_sort,error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="process_inhomogenous",&
         description="If enabled, inhomogenous stacks are also processed by the acc driver.",&
         default_l_val=accdrv_do_inhomogenous,error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="binning_nbins",&
         description="Number of bins used when ordering "&
         //"the stacks with the binning-scheme.",&
         default_i_val=accdrv_binning_nbins,error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="binning_binsize",&
         description="Size of bins used when ordering "&
         //"the stacks with the binning-scheme.",&
         default_i_val=accdrv_binning_binsize,error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL section_add_subsection(section, subsection,error=error)
    CALL section_release(subsection, error=error)
    !
  END SUBROUTINE create_dbcsr_section

! *****************************************************************************
!> \brief   Creates the dbcsr section for configuring FM
!> \param section ...
!> \param error ...
!> \date    2011-04-05
!> \author  Florian Schiffmann
! *****************************************************************************
  SUBROUTINE create_fm_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_fm_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="FM",&
         description="Configuration options for the full matrices.",&
         n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
         error=error)

    NULLIFY (keyword)

    CALL keyword_create(keyword, name="NROW_BLOCKS",&
         description="Defines the number of rows per scalapack block in "//&
         "the creation of block cyclic dense matrices ",&
         default_i_val=32,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="NCOL_BLOCKS",&
         description="Defines the number of columns per scalapack block in "//&
         "the creation of vlock cyclic dense matrices ",&
         default_i_val=32,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="FORCE_BLOCK_SIZE",&
         description="Ensure for small matrices that the layout is compatible "//&
                     "with bigger ones, i.e. no subdivision is performed (can break LAPACK!!!).",&
         usage="FORCE_BLOCK_SIZE",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="TYPE_OF_MATRIX_MULTIPLICATION",&
         description="Allows to switch between scalapack pdgemm and dbcsr_multiply. "//&
                     "On normal systems pdgemm is recommended on system with GPU "//&
                     "is optimized and can give better performance. NOTE: if DBCSR is employed "//&
                     "FORCE_BLOCK_SIZE should be set. The perfomance on GPU's dependes "//&
                     "crucially on the BLOCK_SIZES. Make sure optimized kernels are available.",&
         usage="TYPE_OF_MATRIX_MULTIPLICATION ELPA",&
         default_i_val=do_pdgemm, &
         enum_i_vals=(/do_pdgemm,do_dbcsr/),&
         enum_c_vals=s2a("PDGEMM","DBCSR_MM"),&
         enum_desc=s2a("Standard scalapack: pdgemm",&
                       "DBCSR_MM is employed. This needs local transformation of the matrices"),&
         supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    !
  END SUBROUTINE create_fm_section


END MODULE input_cp2k_global
