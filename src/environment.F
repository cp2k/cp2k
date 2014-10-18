!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Sets up and terminates the global environment variables
!> \par History
!>      - Merged with Quickstep MODULE start_program_run (17.01.2002,MK)
!>      - Compile information added (16.01.2002,MK)
!>      - Merged with MODULE cp2k_input, some rearrangements (30.10.2002,MK)
!> \author JGH,MK
! *****************************************************************************
MODULE environment
  USE bibliography,                    ONLY: Frigo2005,&
                                             cite_reference
  USE cp2k_info,                       ONLY: &
       compile_arch, compile_date, compile_host, compile_revision, cp2k_home, &
       cp2k_version, cp2k_year, get_runtime_info, id_cp2k_version, &
       r_host_name, r_pid, r_user_name
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_fm_diag,                      ONLY: diag_init
  USE cp_fm_struct,                    ONLY: cp_fm_struct_config
  USE cp_fm_types,                     ONLY: cp_fm_setup
  USE cp_ma_interface,                 ONLY: &
       cp_ma_config, cp_ma_current_thread_run, cp_ma_mempol, cp_ma_mpi_sched, &
       cp_ma_print_machine, cp_ma_print_strategy, cp_ma_run_on, &
       cp_ma_thread_run_on, cp_ma_thread_sched
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr,&
                                             debug_print_level,&
                                             high_print_level,&
                                             low_print_level,&
                                             medium_print_level,&
                                             silent_print_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE fft_tools,                       ONLY: fft3d,&
                                             finalize_fft,&
                                             init_fft
  USE force_env_types,                 ONLY: multiple_fe_list
  USE gamma,                           ONLY: deallocate_md_ftable
  USE global_types,                    ONLY: global_environment_type
  USE header,                          ONLY: cp2k_footer,&
                                             cp2k_header
  USE input_constants,                 ONLY: &
       callgraph_all, callgraph_none, do_cp2k, do_diag_elpa, do_diag_sl, &
       do_diag_sl2, do_eip, do_farming, do_fft_fftw3, do_fft_sg, do_fist, &
       do_qs, do_test, energy_run, id_development_version, mol_dyn_run, &
       none_run
  USE input_cp2k_global,               ONLY: create_global_section
  USE input_enumeration_types,         ONLY: enum_i2c,&
                                             enumeration_type
  USE input_keyword_types,             ONLY: keyword_get,&
                                             keyword_type
  USE input_section_types,             ONLY: &
       section_get_keyword, section_get_lval, section_get_rval, &
       section_release, section_type, section_vals_get, &
       section_vals_get_subs_vals, section_vals_get_subs_vals3, &
       section_vals_type, section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp,&
                                             int_8,&
                                             int_size,&
                                             print_kind_info
  USE machine,                         ONLY: flush_should_flush,&
                                             m_energy,&
                                             m_flush_internal,&
                                             m_memory_details,&
                                             m_procrun
  USE machine_architecture_types,      ONLY: group,&
                                             has_ma,&
                                             local,&
                                             threads
  USE message_passing,                 ONLY: add_mp_perf_env,&
                                             describe_mp_perf_env,&
                                             mp_max,&
                                             mp_sum,&
                                             rm_mp_perf_env
  USE orbital_pointers,                ONLY: deallocate_orbital_pointers
  USE orbital_transformation_matrices, ONLY: deallocate_spherical_harmonics
  USE parallel_rng_types,              ONLY: GAUSSIAN,&
                                             check_rng,&
                                             create_rng_stream,&
                                             init_rng,&
                                             write_rng_matrices,&
                                             write_rng_stream
  USE physcon,                         ONLY: write_physcon
  USE reference_manager,               ONLY: collect_citations_from_ranks,&
                                             print_all_references,&
                                             print_format_journal
  USE string_utilities,                ONLY: ascii_to_string,&
                                             integer_to_string,&
                                             string_to_ascii
  USE termination,                     ONLY: set_error_unit,&
                                             stop_memory,&
                                             stop_program
  USE timings,                         ONLY: add_timer_env,&
                                             rm_timer_env,&
                                             timeset,&
                                             timestop,&
                                             timings_setup_tracing
  USE timings_report,                  ONLY: cost_type_energy,&
                                             cost_type_time,&
                                             timings_report_callgraph,&
                                             timings_report_print

  !$ USE OMP_LIB
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'environment'

  ! *** Public subroutines ***
  PUBLIC :: cp2k_finalize, cp2k_init, cp2k_read, cp2k_setup, cp2k_get_walltime

CONTAINS

! *****************************************************************************
!> \brief Initializes a CP2K run (setting of the global environment variables)
!> \param para_env ...
!> \param output_unit ...
!> \param globenv ...
!> \param input_file_name ...
!> \param wdir ...
!> \par History
!>      JGH (28.11.2001) : default for pp_library_path
!>      - print keys added (17.01.2002, MK)
!>      - merged with cp2k_input (30.10.2002,MK)
!> \author JGH,MK
! *****************************************************************************
  SUBROUTINE cp2k_init(para_env,output_unit,globenv,input_file_name,wdir)

    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: output_unit
    TYPE(global_environment_type), POINTER   :: globenv
    CHARACTER(LEN=*)                         :: input_file_name
    CHARACTER(LEN=*), OPTIONAL :: wdir

    CHARACTER(len=*), PARAMETER :: routineN = 'cp2k_init', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: dev_flag
    INTEGER                                  :: my_output_unit
    TYPE(cp_logger_type), POINTER            :: logger
#if defined(__HPM)
#include "f_hpm.h"
#endif

    ! create a timer_env
    CALL add_timer_env()

    ! message passing performance
    CALL add_mp_perf_env()

    ! Set the default logical output and error unit number
    ! write to screen in case of error on non-io nodes, and hope this is appropriate
    IF (para_env%ionode) THEN
       CALL set_error_unit(output_unit)
    ELSE
       CALL set_error_unit(6)
    ENDIF

    ! Set flag if this is a development version
    dev_flag = ""
    IF (id_cp2k_version==id_development_version) dev_flag=" (Development Version)"

    ! Init the default logger
    IF (para_env%source==para_env%mepos) THEN
       my_output_unit=output_unit
    ELSE
       my_output_unit=-1
    END IF
    NULLIFY(logger)
    CALL cp_logger_create(logger,para_env=para_env,&
         default_global_unit_nr=output_unit, &
         close_global_unit_on_dealloc=.FALSE.)
    CALL cp_add_default_logger(logger)
    CALL cp_logger_release(logger)

    ! Initialize timing
    CALL timeset ("CP2K",globenv%handle)

    ! Print header
    CALL cp2k_header(my_output_unit,wdir)

    IF (my_output_unit>0) THEN
       WRITE (UNIT=my_output_unit,FMT="(/,T2,A,T31,A50)")&
            "CP2K| version string: ",&
            ADJUSTR(TRIM(cp2k_version)//TRIM(dev_flag))
       WRITE (UNIT=my_output_unit,FMT="(T2,A,T41,A40)")&
            "CP2K| source code revision number:",&
            ADJUSTR(compile_revision)
       WRITE (UNIT=my_output_unit,FMT="(T2,A,T41,A40)")&
            "CP2K| is freely available from ",&
            ADJUSTR(TRIM(cp2k_home))
       WRITE (UNIT=my_output_unit,FMT="(T2,A,T31,A50)")&
            "CP2K| Program compiled at",&
            ADJUSTR(compile_date(1:MIN(50,LEN(compile_date))))
       WRITE (UNIT=my_output_unit,FMT="(T2,A,T31,A50)")&
            "CP2K| Program compiled on",&
            ADJUSTR(compile_host(1:MIN(50,LEN(compile_host))))
       WRITE (UNIT=my_output_unit,FMT="(T2,A,T31,A50)")&
            "CP2K| Program compiled for",&
            ADJUSTR(compile_arch(1:MIN(50,LEN(compile_arch))))
       WRITE (UNIT=my_output_unit,FMT="(T2,A,T31,A50)")&
            "CP2K| Input file name",&
            ADJUSTR(TRIM(input_file_name))
       CALL m_flush_internal(my_output_unit)
    END IF

  END SUBROUTINE cp2k_init

! *****************************************************************************
!> \brief echos the list of host names and pids
!> \param para_env ...
!> \param output_unit ...
! *****************************************************************************
  SUBROUTINE echo_all_hosts(para_env,output_unit)
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: output_unit

    CHARACTER(len=*), PARAMETER :: routineN = 'echo_all_hosts', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: string
    INTEGER                                  :: ipe, istat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: all_pid
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: all_host

!   *** Print a list of all started processes ***

    ALLOCATE (all_pid(para_env%num_pe),STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,__LINE__,"all_pid",&
            para_env%num_pe*int_size)
    END IF
    all_pid(:) = 0
    all_pid(para_env%mepos+1) = r_pid

    CALL mp_sum(all_pid,para_env%group)
    ALLOCATE (all_host(30,para_env%num_pe),STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,__LINE__,"all_host",&
            30*para_env%num_pe*int_size)
    END IF
    all_host(:,:) = 0
    CALL string_to_ascii(r_host_name,all_host(:,para_env%mepos+1))
    CALL mp_sum(all_host,para_env%group)
    IF (output_unit>0) THEN

       WRITE (UNIT=output_unit,FMT="(T2,A)") ""
       DO ipe=1,para_env%num_pe
          CALL ascii_to_string(all_host(:,ipe),string)
          WRITE (UNIT=output_unit,FMT="(T2,A,T63,I8,T71,I10)")&
               TRIM(r_user_name)//"@"//TRIM(string)//&
               " has created process number",ipe-1,all_pid(ipe)
       END DO
       WRITE (UNIT=output_unit,FMT="(T2,A)") ""
    END IF
    DEALLOCATE (all_pid,STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,__LINE__,"all_pid")
    END IF
    DEALLOCATE (all_host,STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,__LINE__,"all_host")

    END IF

  END SUBROUTINE echo_all_hosts

! *****************************************************************************
!> \brief echos the list the number of process per host
!> \param para_env ...
!> \param output_unit ...
! *****************************************************************************
  SUBROUTINE echo_all_process_host(para_env,output_unit)
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: output_unit

    CHARACTER(len=*), PARAMETER :: routineN = 'echo_all_process_host', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: string, string_sec
    INTEGER                                  :: ipe, istat, jpe, nr_occu
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: all_pid
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: all_host

    ALLOCATE (all_host(30,para_env%num_pe),STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,__LINE__,"all_host",&
            30*para_env%num_pe*int_size)
    END IF
    all_host(:,:) = 0

     IF (m_procrun(r_pid) .EQ. 1) THEN
       CALL string_to_ascii(r_host_name,all_host(:,para_env%mepos+1))
       CALL mp_sum(all_host,para_env%group)
     ENDIF

    IF (output_unit>0) THEN
      ALLOCATE (all_pid(para_env%num_pe),STAT=istat)
      IF (istat /= 0) THEN
         CALL stop_memory(routineN,moduleN,__LINE__,"all_pid",&
            para_env%num_pe*int_size)
      END IF
      all_pid(:) = 0

      WRITE (UNIT=output_unit,FMT="(T2,A)") ""
       DO ipe=1,para_env%num_pe
          nr_occu = 0
        IF (all_pid(ipe) .NE. -1) THEN
          CALL ascii_to_string(all_host(:,ipe),string)
          DO jpe=1,para_env%num_pe
             CALL ascii_to_string(all_host(:,jpe),string_sec)
             IF  (string .EQ. string_sec) THEN
                nr_occu = nr_occu +1
                all_pid(jpe) = -1
             ENDIF
          END DO
          WRITE (UNIT=output_unit,FMT="(T2,A,T63,I8,A)")&
               TRIM(r_user_name)//"@"//TRIM(string)//&
               " is running ",nr_occu, " processes"
          WRITE (UNIT=output_unit,FMT="(T2,A)") ""
        END IF
       END DO
       DEALLOCATE (all_pid,STAT=istat)
       IF (istat /= 0) THEN
         CALL stop_memory(routineN,moduleN,__LINE__,"all_pid")
       END IF

    END IF

    DEALLOCATE (all_host,STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,__LINE__,"all_host")
    END IF

  END SUBROUTINE echo_all_process_host



! *****************************************************************************
!> \brief read part of cp2k_init
!> \param root_section ...
!> \param para_env ...
!> \param globenv the globenv
!> \param error ...
!> \author fawzi
!> \note
!>      The following routines need to be synchronized wrt. adding/removing
!>      of the default environments (logging, perormance,error):
!>      environment:cp2k_init, environment:cp2k_finalize,
!>      f77_interface:f_env_add_defaults, f77_interface:f_env_rm_defaults,
!>      f77_interface:create_force_env, f77_interface:destroy_force_env
! *****************************************************************************
  SUBROUTINE cp2k_read(root_section,para_env,globenv,error)

    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp2k_read', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: c_val
    INTEGER                                  :: iw
    TYPE(cp_logger_type), POINTER            :: logger

!   *** Read the input/output section ***

    logger => cp_get_default_logger()
#if defined(__HPM)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         root_section,"GLOBAL%PRINT/HPM",error=error),cp_p_file)) THEN
       CALL f_hpminit(para_env%mepos,"CP2K")
       use_HPM = .TRUE.
    ENDIF
#endif
    ! try to use better names for the local log if it is not too late
    CALL section_vals_val_get(root_section,"GLOBAL%OUTPUT_FILE_NAME",&
         c_val=c_val,error=error)
    IF (c_val/="") THEN
       CALL cp_logger_set(logger,&
            local_filename=TRIM(c_val)//"_localLog")
    END IF
    CALL section_vals_val_get(root_section,"GLOBAL%PROJECT",c_val=c_val,error=error)
    IF (c_val/="") THEN
       CALL cp_logger_set(logger,local_filename=TRIM(c_val)//"_localLog")
    END IF
    logger%iter_info%project_name=c_val
    CALL section_vals_val_get(root_section,"GLOBAL%PRINT_LEVEL",i_val=logger%iter_info%print_level,error=error)

    !   *** Read the CP2K section ***
    CALL read_cp2k_section(root_section,para_env,globenv,error=error)

    iw = cp_print_key_unit_nr(logger,root_section,"GLOBAL%PRINT/BASIC_DATA_TYPES",&
         extension=".Log",error=error)
    IF (iw>0) CALL print_kind_info(iw)
    CALL cp_print_key_finished_output(iw,logger,root_section,&
         "GLOBAL%PRINT/BASIC_DATA_TYPES",error=error)

    iw = cp_print_key_unit_nr(logger,root_section,"GLOBAL%PRINT/PHYSCON",&
         extension=".Log",error=error)
    IF (iw>0) CALL write_physcon(iw)
    CALL cp_print_key_finished_output(iw,logger,root_section,&
         "GLOBAL%PRINT/PHYSCON",error=error)

  END SUBROUTINE cp2k_read

! *****************************************************************************
!> \brief globenv initializations that need the input and error
!> \param root_section ...
!> \param para_env ...
!> \param globenv the global environment to initialize
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
!> \note
!>      if possible do the initializations here as the environement
!>      (error,...) is setup, instaed of cp2k_init
! *****************************************************************************
  SUBROUTINE cp2k_setup(root_section,para_env,globenv,error)

    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp2k_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: input_seed, iw
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(3, 2)           :: initial_seed
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    IF (.NOT.failure) THEN

       ! Initialize the parallel random number generator

       CALL init_rng()
       iw = cp_print_key_unit_nr(logger,root_section,"GLOBAL%PRINT/RNG_MATRICES",&
            extension=".Log",error=error)
       IF (iw > 0) THEN
          CALL write_rng_matrices(iw)
       END IF
       CALL cp_print_key_finished_output(iw,logger,root_section,&
            "GLOBAL%PRINT/RNG_MATRICES",&
            error=error)

       ! Initialize a global normally Gaussian distributed (pseudo)random number stream

       CALL section_vals_val_get(root_section,"GLOBAL%SEED",i_val=input_seed,error=error)
       initial_seed(:,:) = REAL(input_seed,KIND=dp)
       CALL create_rng_stream(rng_stream=globenv%gaussian_rng_stream,&
            name="Global Gaussian random numbers",&
            distribution_type=GAUSSIAN,&
            seed=initial_seed,&
            extended_precision=.TRUE.,&
            error=error)

       iw = cp_print_key_unit_nr(logger,root_section,"GLOBAL%PRINT/RNG_CHECK",&
            extension=".Log",error=error)
       IF (iw > 0) THEN
          CALL check_rng(iw,para_env%ionode,error)
       END IF
       CALL cp_print_key_finished_output(iw,logger,root_section,&
            "GLOBAL%PRINT/RNG_CHECK",&
            error=error)

       iw = cp_print_key_unit_nr(logger,root_section,"GLOBAL%PRINT/GLOBAL_GAUSSIAN_RNG",&
            extension=".Log",error=error)
       IF (iw > 0) THEN
          CALL write_rng_stream(globenv%gaussian_rng_stream,iw,write_all=.TRUE.,error=error)
       END IF
       CALL cp_print_key_finished_output(iw,logger,root_section,&
            "GLOBAL%PRINT/GLOBAL_GAUSSIAN_RNG",&
            error=error)

    END IF

  END SUBROUTINE cp2k_setup

! *****************************************************************************
!> \brief read the global sectionof new input
!> \param root_section ...
!> \param para_env ...
!> \param globenv ...
!> \param error ...
!> \par History
!>      06-2005 [created]
!> \author MI
!> \note
!>      Should not be required anymore once everything is converted
!>      to get information directly from the input structure
! *****************************************************************************
  SUBROUTINE read_global_section(root_section,para_env,globenv,error)

    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_global_section', &
      routineP = moduleN//':'//routineN, start_section_label = "GLOBAL"

    CHARACTER(len=13)                        :: tracing_string
    CHARACTER(len=6)                         :: print_level_string
    CHARACTER(len=default_path_length) :: basis_set_file_name, &
      coord_file_name, geminal_file_name, mm_potential_file_name, &
      potential_file_name
    CHARACTER(len=default_string_length)     :: env_num, project_name
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: trace_routines
    INTEGER :: i_diag, i_fft, iforce_eval, method_name_id, n_rep_val, &
      nforce_eval, num_threads, output_unit, print_level, trace_max, unit_nr

!$  INTEGER :: nid
    INTEGER(kind=int_8) :: Buffers, Buffers_avr, Buffers_max, Buffers_min, &
      Cached, Cached_avr, Cached_max, Cached_min, MemFree, MemFree_avr, &
      MemFree_max, MemFree_min, MemLikelyFree, MemLikelyFree_avr, &
      MemLikelyFree_max, MemLikelyFree_min, MemTotal, MemTotal_avr, &
      MemTotal_max, MemTotal_min, Slab, Slab_avr, Slab_max, Slab_min, &
      SReclaimable, SReclaimable_avr, SReclaimable_max, SReclaimable_min
    INTEGER, DIMENSION(:), POINTER           :: i_force_eval
    LOGICAL                                  :: ata, do_echo_all_hosts, efl, &
                                                explicit, failure, trace, &
                                                trace_master
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(enumeration_type), POINTER          :: enum1, enum2
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section
    TYPE(section_vals_type), POINTER         :: dft_section, &
                                                force_env_sections, &
                                                global_section, qmmm_section, &
                                                subsys_section

    failure = .FALSE.
    NULLIFY(dft_section,global_section, i_force_eval)

    global_section => section_vals_get_subs_vals(root_section,"GLOBAL",error=error)
    CALL section_vals_val_get(global_section,"BLACS_GRID",i_val=globenv%blacs_grid_layout,error=error)
    CALL section_vals_val_get(global_section,"BLACS_REPEATABLE",l_val=globenv%blacs_repeatable,error=error)
    CALL section_vals_val_get(global_section,"PREFERRED_DIAG_LIBRARY",i_val=i_diag,error=error)
    CALL section_vals_val_get(global_section,"PREFERRED_FFT_LIBRARY",i_val=i_fft,error=error)

    CALL section_vals_val_get(global_section,"PRINT_LEVEL",i_val=print_level,error=error)
    CALL section_vals_val_get(global_section,"PROGRAM_NAME",i_val=globenv%prog_name_id,error=error)
    CALL section_vals_val_get(global_section,"FFT_POOL_SCRATCH_LIMIT",i_val=globenv%fft_pool_scratch_limit,error=error)
    CALL section_vals_val_get(global_section,"FFTW_PLAN_TYPE",i_val=globenv%fftw_plan_type,error=error)
    CALL section_vals_val_get(global_section,"PROJECT_NAME",c_val=project_name,error=error)
    CALL section_vals_val_get(global_section,"FFTW_WISDOM_FILE_NAME",c_val=globenv%fftw_wisdom_file_name,error=error)
    CALL section_vals_val_get(global_section,"RUN_TYPE",i_val=globenv%run_type_id,error=error)
    CALL cp2k_get_walltime(section=global_section, keyword_name="WALLTIME",&
                           walltime=globenv%cp2k_target_time, error=error)
    CALL section_vals_val_get(global_section,"TRACE",l_val=trace,error=error)
    CALL section_vals_val_get(global_section,"TRACE_MASTER",l_val=trace_MASTER,error=error)
    CALL section_vals_val_get(global_section,"TRACE_MAX",i_val=trace_max,error=error)
    CALL section_vals_val_get(global_section,"TRACE_ROUTINES",explicit=explicit,error=error)
    IF (explicit) THEN
       CALL section_vals_val_get(global_section,"TRACE_ROUTINES",c_vals=trace_routines,error=error)
    ELSE
       NULLIFY(trace_routines)
    ENDIF
    CALL section_vals_val_get(global_section,"FLUSH_SHOULD_FLUSH",l_val=flush_should_flush,error=error)
    CALL section_vals_val_get(global_section,"ECHO_ALL_HOSTS",l_val=do_echo_all_hosts,error=error)
    force_env_sections => section_vals_get_subs_vals(root_section,"FORCE_EVAL",error=error)
    CALL section_vals_get(force_env_sections, n_repetition=nforce_eval, error=error)
    logger=>cp_error_get_logger(error)
    output_unit=cp_print_key_unit_nr(logger,global_section,"PROGRAM_RUN_INFO",&
         extension=".log",error=error)
    
    CALL fm_setup(global_section,error)
     
    IF(trace .AND. (.NOT. trace_master .OR. para_env%mepos==0)) THEN
       unit_nr = -1
       IF (logger%para_env%mepos==logger%para_env%source .OR. .NOT. trace_master) &
           unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
       WRITE (tracing_string,"(I6.6,A1,I6.6)") para_env%mepos,":",para_env%num_pe
       IF(ASSOCIATED(trace_routines))  THEN
          CALL timings_setup_tracing(trace_max, unit_nr, tracing_string, trace_routines)
       ELSE
          CALL timings_setup_tracing(trace_max, unit_nr, tracing_string)
       END IF
    ENDIF 

    SELECT CASE(i_diag)
    CASE(do_diag_sl)
       globenv%diag_library="SL" 
    CASE(do_diag_sl2)
       globenv%diag_library="SL2" 
    CASE(do_diag_elpa)
       globenv%diag_library="ELPA" 
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,&
                         "Unknown DIAG type",para_env)
    END SELECT

    SELECT CASE(i_fft)
    CASE(do_fft_sg)
       globenv%default_fft_library="FFTSG"
    CASE(do_fft_fftw3)
       globenv%default_fft_library="FFTW3"
       CALL cite_reference(Frigo2005)
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,&
                         "Unknown FFT type",para_env)
    END SELECT

    IF (globenv%run_type_id == 0) THEN
       SELECT CASE (globenv%prog_name_id)
       CASE (do_farming, do_test)
          globenv%run_type_id = none_run
       CASE (do_cp2k)
          IF (nforce_eval /=1) THEN
             ! multiple force_eval corresponds at the moment to RESPA calculations only
             ! default MD
             globenv%run_type_id = mol_dyn_run
          ELSE
             CALL section_vals_val_get(force_env_sections,"METHOD",i_val=method_name_id,error=error)
             SELECT CASE (method_name_id)
             CASE (do_fist)
                globenv%run_type_id = mol_dyn_run
             CASE (do_eip)
                globenv%run_type_id = mol_dyn_run
             CASE (do_qs)
                globenv%run_type_id = energy_run
             END SELECT
          END IF
       END SELECT
    END IF

    IF(globenv%prog_name_id == do_farming .AND. globenv%run_type_id /= none_run) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
                         "FARMING program supports only NONE as run type",para_env)
    ENDIF

    IF (globenv%prog_name_id == do_test .AND. globenv%run_type_id /= none_run) &
       CALL stop_program(routineN,moduleN,__LINE__,&
                         "TEST program supports only NONE as run type",para_env)

    CALL m_memory_details(MemTotal,MemFree,Buffers,Cached,Slab,SReclaimable,MemLikelyFree)
    MemTotal_avr=MemTotal
    MemFree_avr=MemFree
    Buffers_avr=Buffers
    Cached_avr=Cached
    Slab_avr=Slab
    SReclaimable_avr=SReclaimable
    MemLikelyFree_avr=MemLikelyFree
    CALL mp_sum(MemTotal_avr,para_env%group) ; MemTotal_avr=MemTotal_avr/para_env%num_pe/1024
    CALL mp_sum(MemFree_avr,para_env%group) ; MemFree_avr=MemFree_avr/para_env%num_pe/1024
    CALL mp_sum(Buffers_avr,para_env%group) ; Buffers_avr=Buffers_avr/para_env%num_pe/1024
    CALL mp_sum(Cached_avr,para_env%group) ; Cached_avr=Cached_avr/para_env%num_pe/1024
    CALL mp_sum(Slab_avr,para_env%group) ; Slab_avr=Slab_avr/para_env%num_pe/1024
    CALL mp_sum(SReclaimable_avr,para_env%group) ; SReclaimable_avr=SReclaimable_avr/para_env%num_pe/1024
    CALL mp_sum(MemLikelyFree_avr,para_env%group) ; MemLikelyFree_avr=MemLikelyFree_avr/para_env%num_pe/1024

    MemTotal_min=-MemTotal
    MemFree_min=-MemFree
    Buffers_min=-Buffers
    Cached_min=-Cached
    Slab_min=-Slab
    SReclaimable_min=-SReclaimable
    MemLikelyFree_min=-MemLikelyFree
    CALL mp_max(MemTotal_min,para_env%group) ; MemTotal_min=-MemTotal_min/1024
    CALL mp_max(MemFree_min,para_env%group) ; MemFree_min=-MemFree_min/1024
    CALL mp_max(Buffers_min,para_env%group) ; Buffers_min=-Buffers_min/1024
    CALL mp_max(Cached_min,para_env%group) ; Cached_min=-Cached_min/1024
    CALL mp_max(Slab_min,para_env%group) ; Slab_min=-Slab_min/1024
    CALL mp_max(SReclaimable_min,para_env%group) ; SReclaimable_min=-SReclaimable_min/1024
    CALL mp_max(MemLikelyFree_min,para_env%group) ; MemLikelyFree_min=-MemLikelyFree_min/1024

    MemTotal_max=MemTotal
    MemFree_max=MemFree
    Buffers_max=Buffers
    Cached_max=Cached
    Slab_max=Slab
    SReclaimable_max=SReclaimable
    MemLikelyFree_max=MemLikelyFree
    CALL mp_max(MemTotal_max,para_env%group) ; MemTotal_max=MemTotal_max/1024
    CALL mp_max(MemFree_max,para_env%group) ; MemFree_max=MemFree_max/1024
    CALL mp_max(Buffers_max,para_env%group) ; Buffers_max=Buffers_max/1024
    CALL mp_max(Cached_max,para_env%group) ; Cached_max=Cached_max/1024
    CALL mp_max(Slab_max,para_env%group) ; Slab_max=Slab_max/1024
    CALL mp_max(SReclaimable_max,para_env%group) ; SReclaimable_max=SReclaimable_max/1024
    CALL mp_max(MemLikelyFree_max,para_env%group) ; MemLikelyFree_max=MemLikelyFree_max/1024

    MemTotal=MemTotal/1024
    MemFree=MemFree/1024
    Buffers=Buffers/1024
    Cached=Cached/1024
    Slab=Slab/1024
    SReclaimable=SReclaimable/1024
    MemLikelyFree=MemLikelyFree/1024



   ! Initialize the MA configurations
    IF (has_ma) THEN
      CALL cp_ma_config(root_section,error)
      CALL cp_ma_mempol(error)
    ENDIF

    !   *** Print a list of all started processes ***
    IF (do_echo_all_hosts) THEN
     CALL echo_all_hosts(para_env,output_unit)

     ! *** Print the number of processes per host ***
     CALL echo_all_process_host(para_env,output_unit)
    ENDIF

    num_threads=1
!$  num_threads=omp_get_max_threads()
    IF (output_unit > 0 ) THEN
       WRITE (UNIT=output_unit,FMT=*)
       CALL multiple_fe_list(force_env_sections, root_section, i_force_eval, nforce_eval, error)
       DO iforce_eval = 1, nforce_eval
          dft_section => section_vals_get_subs_vals3(force_env_sections,"DFT",&
               i_rep_section=i_force_eval(iforce_eval),error=error)
          qmmm_section => section_vals_get_subs_vals3(force_env_sections,"QMMM",&
               i_rep_section=i_force_eval(iforce_eval),error=error)
          CALL section_vals_val_get(dft_section,"BASIS_SET_FILE_NAME",&
               c_val=basis_set_file_name , ignore_required=.TRUE., error=error)
          CALL section_vals_val_get(dft_section,"GEMINAL_FILE_NAME",&
               c_val=geminal_file_name , ignore_required=.TRUE., error=error)
          CALL section_vals_val_get(dft_section,"POTENTIAL_FILE_NAME",&
               c_val=potential_file_name , ignore_required=.TRUE., error=error)

          CALL section_vals_val_get(qmmm_section,"MM_POTENTIAL_FILE_NAME",&
               c_val=mm_potential_file_name, ignore_required=.TRUE., error=error)
          ! SUBSYS - If any
          subsys_section => section_vals_get_subs_vals3(force_env_sections,"SUBSYS",&
               i_rep_section=i_force_eval(iforce_eval),error=error)
          CALL section_vals_get(subsys_section, explicit=explicit, error=error)
          coord_file_name = "__STD_INPUT__"
          IF (explicit) THEN
             CALL section_vals_val_get(subsys_section,"TOPOLOGY%COORD_FILE_NAME",&
                  n_rep_val=n_rep_val, error=error)
             IF (n_rep_val==1) THEN
                CALL section_vals_val_get(subsys_section,"TOPOLOGY%COORD_FILE_NAME",&
                     c_val=coord_file_name, error=error)
             END IF
          END IF
          CALL integer_to_string(i_force_eval(iforce_eval),env_num)

          WRITE (UNIT=output_unit,FMT="(T2,A,T41,A)")&
               start_section_label//"| Force Environment number",&
               ADJUSTR(env_num(:40)),&
               start_section_label//"| Basis set file name",&
               ADJUSTR(basis_set_file_name(:40)),&
               start_section_label//"| Geminal file name",&
               ADJUSTR(geminal_file_name(:40)),&
               start_section_label//"| Potential file name",&
               ADJUSTR(potential_file_name(:40)),&
               start_section_label//"| MM Potential file name",&
               ADJUSTR(mm_potential_file_name(:40)),&
               start_section_label//"| Coordinate file name",&
               ADJUSTR(coord_file_name(:40))
       END DO
       DEALLOCATE(i_force_eval)

       NULLIFY (enum1,enum2,keyword,section)
       CALL create_global_section(section,error=error)
       keyword => section_get_keyword(section,"PROGRAM_NAME",error=error)
       CALL keyword_get(keyword,enum=enum1,error=error)
       keyword => section_get_keyword(section,"RUN_TYPE",error=error)
       CALL keyword_get(keyword,enum=enum2,error=error)

       WRITE (UNIT=output_unit,FMT="(T2,A,T41,A40)")&
            start_section_label//"| Method name",&
            ADJUSTR(TRIM(enum_i2c(enum1,globenv%prog_name_id,error=error))),&
            start_section_label//"| Project name",&
            ADJUSTR(project_name(:40)),&
            start_section_label//"| Preferred FFT library",&
            ADJUSTR(globenv%default_fft_library(:40)),&
            start_section_label//"| Preferred diagonalization lib.",&
            ADJUSTR(globenv%diag_library(:40)),&
            start_section_label//"| Run type",&
            ADJUSTR(TRIM(enum_i2c(enum2,globenv%run_type_id,error=error)))

       CALL section_release(section,error=error)

       CALL section_vals_val_get(global_section,"ALLTOALL_SGL",l_val=ata,error=error)
       WRITE (UNIT=output_unit,FMT="(T2,A,T80,L1)")&
            start_section_label//"| All-to-all communication in single precision",ata
       CALL section_vals_val_get(global_section,"EXTENDED_FFT_LENGTHS",l_val=efl,error=error)
       WRITE (UNIT=output_unit,FMT="(T2,A,T80,L1)")&
            start_section_label//"| FFTs using library dependent lengths",efl

       SELECT CASE(print_level)
       CASE(silent_print_level)
         print_level_string="SILENT"
       CASE(low_print_level)
         print_level_string="   LOW"
       CASE(medium_print_level)
         print_level_string="MEDIUM"
       CASE(high_print_level)
         print_level_string="  HIGH"
       CASE(debug_print_level)
         print_level_string=" DEBUG"
       CASE DEFAULT
         CALL stop_program(routineN,moduleN,__LINE__,&
                           "Unknown print_level",para_env)
       END SELECT

       WRITE (UNIT=output_unit,FMT="(T2,A,T75,A6)")&
            start_section_label//"| Global print level",print_level_string
       WRITE (UNIT=output_unit,FMT="(T2,A,T75,I6)")&
            start_section_label//"| Total number of message passing processes",&
            para_env%num_pe,&
            start_section_label//"| Number of threads for this process",&
            num_threads,&
            start_section_label//"| This output is from process",para_env%mepos

       WRITE (UNIT=output_unit,FMT="()")
       WRITE (UNIT=output_unit,FMT="(T2,A)") "MEMORY| system memory details [Kb]"
       WRITE (UNIT=output_unit,FMT="(T2,A23,4A14)") "MEMORY|                ","rank 0", "min", "max", "average"
       WRITE (UNIT=output_unit,FMT="(T2,A23,4I14)") "MEMORY| MemTotal       ",memtotal,memtotal_min,memtotal_max,memtotal_avr
       WRITE (UNIT=output_unit,FMT="(T2,A23,4I14)") "MEMORY| MemFree        ",memFree,memfree_min,memfree_max,memfree_avr
       WRITE (UNIT=output_unit,FMT="(T2,A23,4I14)") "MEMORY| Buffers        ",Buffers,Buffers_min,Buffers_max,Buffers_avr
       WRITE (UNIT=output_unit,FMT="(T2,A23,4I14)") "MEMORY| Cached         ",Cached,Cached_min,Cached_max,Cached_avr
       WRITE (UNIT=output_unit,FMT="(T2,A23,4I14)") "MEMORY| Slab           ",Slab,Slab_min,Slab_max,Slab_avr
       WRITE (UNIT=output_unit,FMT="(T2,A23,4I14)") "MEMORY| SReclaimable   ",SReclaimable,SReclaimable_min,SReclaimable_max,&
                                                             SReclaimable_avr
       WRITE (UNIT=output_unit,FMT="(T2,A23,4I14)") "MEMORY| MemLikelyFree  ",MemLikelyFree,MemLikelyFree_min,MemLikelyFree_max,&
                                                             MemLikelyFree_avr
       WRITE (UNIT=output_unit,FMT='()')

#if defined __GEMINI || __SEASTAR || __BLUEGENE || __NET || __HWLOC || __LIBNUMA
    WRITE(UNIT=output_unit,FMT='(T2,A)') "MACHINE HARDWARE TOPOLOGY| Strategies details"
    CALL cp_ma_print_strategy(unit_nr=output_unit,error=error)
    WRITE (UNIT=output_unit,FMT='()')
#endif
       ! Show machine topology
       CALL cp_ma_print_machine(unit_nr=output_unit, error=error)

    END IF

    IF (has_ma) THEN
      CALL cp_ma_mpi_sched()
      CALL cp_ma_run_on(para_env, error, output_unit)


!$omp parallel private(nid)
!$omp critical
!$          nid = omp_get_thread_num()
!$          CALL cp_ma_thread_sched(nid, error=error )
!$omp end critical
!$omp end parallel

       CALL cp_ma_thread_run_on(output_unit, error=error)

!$omp parallel private(nid)
!$          nid = omp_get_thread_num()
!$          CALL cp_ma_current_thread_run(nid, error=error)
!$omp end parallel
     ENDIF

    CALL cp_print_key_finished_output(output_unit,logger,global_section,&
         "PROGRAM_RUN_INFO", error=error)

  END SUBROUTINE read_global_section

! *****************************************************************************
!> \brief ...
!> \param root_section ...
!> \param para_env ...
!> \param globenv ...
!> \param error ...
!> \par History
!>      2-Dec-2000 (JGH) added default fft library
!> \author JGH,MK
! *****************************************************************************
  SUBROUTINE read_cp2k_section(root_section,para_env,globenv,error)

    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_cp2k_section', &
      routineP = moduleN//':'//routineN, start_section_label = "CP2K"

    INTEGER                                  :: output_unit
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: global_section

    global_section => section_vals_get_subs_vals(root_section,"GLOBAL",error=error)
    CALL read_global_section(root_section,para_env,globenv,error=error)
    logger => cp_error_get_logger(error)
    output_unit=cp_print_key_unit_nr(logger,global_section,"PROGRAM_RUN_INFO",&
         extension=".log",error=error)

    CALL fft_setup_library(globenv, para_env, global_section, output_unit,error=error)
    CALL diag_setup_library(globenv, output_unit,error=error)

    CALL cp_print_key_finished_output(output_unit,logger,global_section,&
         "PROGRAM_RUN_INFO", error=error)

  END SUBROUTINE read_cp2k_section

! *****************************************************************************
!> \brief check FFT preferred library availability, if not switch
!> \param globenv ...
!> \param para_env ...
!> \param global_section ...
!> \param output_unit ...
!> \param error ...
!> \par History
!>      2-Dec-2000 (JGH) added default fft library
!>      Nov-2013 (MI) refactoring
!> \author JGH,MK
! *****************************************************************************


  SUBROUTINE fft_setup_library(globenv, para_env,global_section, output_unit,error)

    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: global_section
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fft_setup_library', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(4, 4, 4)     :: zz
    INTEGER                                  :: stat
    INTEGER, DIMENSION(3)                    :: n
    LOGICAL                                  :: try_fftw

    n(:) = 4
    zz(:,:,:) = 0.0_dp

    !
    ! Setup the FFT library
    ! If the user has specified PREFERRED_FFT_LIBRARY try that first (default FFTW3)
    ! If that one is not available, try FFTW3 (unless it has been tried already)
    ! If FFTW3 is not available use FFTSG
    !
    IF (globenv%default_fft_library.EQ."FFTW3") THEN
       try_fftw=.FALSE.
    ELSE
       try_fftw=.TRUE.
    ENDIF

    !   *** Initialize FFT library with the user's prefered FFT library ***
    CALL init_fft(fftlib=TRIM(globenv%default_fft_library),&
         alltoall=section_get_lval(global_section,"ALLTOALL_SGL",error), &
         fftsg_sizes=.NOT.section_get_lval(global_section,"EXTENDED_FFT_LENGTHS",error), &
         pool_limit=globenv%fft_pool_scratch_limit,&
         wisdom_file=globenv%fftw_wisdom_file_name,&
         plan_style=globenv%fftw_plan_type,&
         error=error)

    !   *** Check for FFT library ***
    CALL fft3d(1,n,zz,status=stat)
    IF (stat /= 0) THEN

       IF (try_fftw) THEN
          IF (output_unit > 0) THEN
             WRITE(output_unit,'(A,A,T55,A)') &
                  " WARNING : FFT library "//TRIM(globenv%default_fft_library)//&
                  " is not available "," Trying FFTW3"
          ENDIF
          globenv%default_fft_library="FFTW3"
          CALL init_fft(fftlib=TRIM(globenv%default_fft_library),&
               alltoall=section_get_lval(global_section,"ALLTOALL_SGL",error), &
               fftsg_sizes=.NOT.section_get_lval(global_section,"EXTENDED_FFT_LENGTHS",error), &
               pool_limit=globenv%fft_pool_scratch_limit,&
               wisdom_file=globenv%fftw_wisdom_file_name,&
               plan_style=globenv%fftw_plan_type,&
               error=error)

          CALL fft3d(1,n,zz,status=stat)
       ENDIF

       IF (stat /=0 ) THEN

          IF (output_unit > 0) THEN
             WRITE(output_unit,'(A,A,T55,A)') &
                  " WARNING : FFT library "//TRIM(globenv%default_fft_library)//&
                  " is not available "," Trying FFTSG as a default"
          ENDIF

          globenv%default_fft_library="FFTSG"
          CALL init_fft(fftlib=TRIM(globenv%default_fft_library),&
               alltoall=section_get_lval(global_section,"ALLTOALL_SGL",error), &
               fftsg_sizes=.NOT.section_get_lval(global_section,"EXTENDED_FFT_LENGTHS",error), &
               pool_limit=globenv%fft_pool_scratch_limit,&
               wisdom_file=globenv%fftw_wisdom_file_name,&
               plan_style=globenv%fftw_plan_type,&
               error=error)

          CALL fft3d(1,n,zz,status=stat)
          IF (stat /= 0) THEN
             CALL stop_program(routineN,moduleN,__LINE__,&
                               "FFTSG not functional....",para_env)
          ENDIF

       ENDIF

    END IF

  END SUBROUTINE fft_setup_library

! *****************************************************************************
!> \brief availability diagonalizatioon library
!>
!> \param globenv ...
!> \param output_unit ...
!> \param error ...
!> \author MI
! *****************************************************************************
  SUBROUTINE diag_setup_library(globenv, output_unit,error)
    TYPE(global_environment_type), POINTER   :: globenv
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'diag_setup_library', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: switched

    switched = .FALSE.
  
    CALL diag_init(diag_lib=TRIM(globenv%diag_library), switched=switched,error=error)
  
    IF(switched) THEN
  
          IF (output_unit > 0) THEN
             WRITE(output_unit,'(A,A,T55,A)') &
                  " WARNING : DIAGONALIZATION library "//TRIM(globenv%diag_library)//&
                  " is not available "," Trying SCALAPACK"
          ENDIF
  
    END IF
  
  END SUBROUTINE diag_setup_library

! *****************************************************************************
!> \brief ...
!> \param glob_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE fm_setup(glob_section,error)
    TYPE(section_vals_type), POINTER         :: glob_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fm_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: multiplication_type, ncb, nrb
    LOGICAL                                  :: force_me
    TYPE(section_vals_type), POINTER         :: fm_section

    fm_section => section_vals_get_subs_vals(glob_section,"FM",error=error)

    CALL section_vals_val_get(fm_section,"NROW_BLOCKS",i_val=nrb,error=error)
    CALL section_vals_val_get(fm_section,"NCOL_BLOCKS",i_val=ncb,error=error)
    CALL section_vals_val_get(fm_section,"FORCE_BLOCK_SIZE",l_val=force_me,error=error)

    CALL cp_fm_struct_config(nrow_block=nrb,ncol_block=ncb,force_block=force_me)

    CALL section_vals_val_get(fm_section,"TYPE_OF_MATRIX_MULTIPLICATION",&
                              i_val=multiplication_type,error=error)

    CALL cp_fm_setup(multiplication_type)

  END SUBROUTINE fm_setup
! *****************************************************************************
!> \brief reads the Walltime also in format HH:MM:SS
!> \param section ...
!> \param keyword_name ...
!> \param walltime ...
!> \param error ...
!> \par History
!>      none
!> \author Mandes
! *****************************************************************************
  SUBROUTINE cp2k_get_walltime(section, keyword_name, walltime, error)
    TYPE(section_vals_type), POINTER         :: section
    CHARACTER(len=*), INTENT(in)             :: keyword_name
    REAL(KIND=dp), INTENT(out)               :: walltime
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp2k_get_walltime', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=100)                       :: ctmp
    INTEGER                                  :: ierr, itmp
    LOGICAL                                  :: failure

    CALL section_vals_val_get(section, keyword_name, c_val=ctmp, error=error)

    IF(ctmp == "") THEN
      walltime = -1.0_dp
    ELSE
      IF(INDEX(ctmp,":",BACK=.TRUE.).GT.0)THEN
        walltime = 0
        ! read seconds
        READ(ctmp(INDEX(ctmp,":",BACK=.TRUE.)+1:),FMT=*,IOSTAT=ierr) walltime
        CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
        ! read minutes
        READ(ctmp(INDEX(ctmp,":",BACK=.TRUE.)-2:INDEX(ctmp,":",BACK=.TRUE.)-1),FMT=*,IOSTAT=ierr) itmp
        walltime = walltime + REAL(60*itmp,KIND=dp)
        ! read hours
        READ(ctmp(1:INDEX(ctmp,":")-1),FMT=*,IOSTAT=ierr) itmp
        CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
        walltime = walltime + REAL(3600*itmp,KIND=dp)
      ELSE
        READ(ctmp,FMT=*,IOSTAT=ierr) walltime
        CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
      END IF
    END IF
  END SUBROUTINE cp2k_get_walltime

! *****************************************************************************
!> \brief Writes final timings and banner for CP2K
!> \param root_section ...
!> \param para_env ...
!> \param globenv ...
!> \param wdir ...
!> \param error ...
!> \par History
!>      none
!> \author JGH,MK
!> \note
!>      The following routines need to be synchronized wrt. adding/removing
!>      of the default environments (logging, perormance,error):
!>      environment:cp2k_init, environment:cp2k_finalize,
!>      f77_interface:f_env_add_defaults, f77_interface:f_env_rm_defaults,
!>      f77_interface:create_force_env, f77_interface:destroy_force_env
! *****************************************************************************
  SUBROUTINE cp2k_finalize(root_section,para_env,globenv,wdir,error)

    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    CHARACTER(LEN=*), OPTIONAL               :: wdir
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp2k_finalize', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: dev_flag
    CHARACTER(LEN=default_path_length)       :: cg_filename
    INTEGER                                  :: iw, unit_exit, cg_mode
    LOGICAL                                  :: delete_it,failure,&
                                                sort_by_self_time
    REAL(KIND=dp)                            :: r_timings
    TYPE(cp_logger_type), POINTER            :: logger

#if defined(__HPM)
#include "f_hpm.h"
#endif
! look if we inherited a failure, more care is needed if so
! i.e. the input is most likely not available

    failure = .FALSE.
    CALL cp_error_check(error,failure)

    ! Set flag if this is a development version
    dev_flag = ""
    IF (id_cp2k_version==id_development_version) dev_flag=" (Development Version)"

    ! Clean up
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    CALL deallocate_spherical_harmonics()
    CALL deallocate_orbital_pointers()
    CALL deallocate_md_ftable()

    ! finalize the fft (i.e. writes the wisdom if FFTW3 )
    CALL finalize_fft(para_env,globenv%fftw_wisdom_file_name, error)

    ! Write message passing performance info

    IF (.NOT. failure) THEN
       iw=cp_print_key_unit_nr(logger,root_section,"GLOBAL%PROGRAM_RUN_INFO",&
            extension=".log",error=error)
       CALL describe_mp_perf_env ( iw )
       CALL cp_print_key_finished_output(iw,logger,root_section,&
            "GLOBAL%PROGRAM_RUN_INFO", error=error)
    ENDIF

#if defined(__HPM)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         root_section,"GLOBAL%PRINT/HPM",error=error),cp_p_file)) THEN
       CALL f_hpmterminate(para_env%group)
       use_HPM=.FALSE.
    ENDIF
#endif

    IF (.NOT. failure) THEN
       CALL collect_citations_from_ranks(para_env)
       iw = cp_print_key_unit_nr(logger,root_section,"GLOBAL%REFERENCES",&
            extension=".Log",error=error)
       IF (iw>0) THEN
          WRITE (UNIT=iw,FMT="(/,T2,A)") REPEAT("-",79)
          WRITE (UNIT=iw,FMT="(T2,A,T80,A)") "-","-"
          WRITE (UNIT=iw,FMT="(T2,A,T30,A,T80,A)") "-","R E F E R E N C E S","-"
          WRITE (UNIT=iw,FMT="(T2,A,T80,A)") "-","-"
          WRITE (UNIT=iw,FMT="(T2,A)") REPEAT("-",79)

          WRITE (UNIT=iw,FMT="(T2,A)") ""
          WRITE (UNIT=iw,FMT="(T2,A)") TRIM(cp2k_version)//TRIM(dev_flag)//", the CP2K developers group ("//TRIM(cp2k_year)//")."
          WRITE (UNIT=iw,FMT="(T2,A)") "CP2K is freely available from "//TRIM(cp2k_home)//" ."

          CALL print_all_references(sorted=.TRUE.,cited_only=.TRUE., &
               FORMAT=print_format_journal,unit=iw)
       ENDIF
       CALL cp_print_key_finished_output(iw,logger,root_section,&
            "GLOBAL%REFERENCES",error=error)
    ENDIF

    CALL timestop(globenv%handle) ! corresponding the "CP2K" in cp2k_init
    
    IF (.NOT.failure) THEN
       iw = cp_print_key_unit_nr(logger,root_section,"GLOBAL%TIMINGS",&
            extension=".Log",error=error)
       r_timings = section_get_rval(root_section,"GLOBAL%TIMINGS%THRESHOLD",error)
       sort_by_self_time = section_get_lval(root_section,"GLOBAL%TIMINGS%SORT_BY_SELF_TIME",error)
       IF (m_energy().NE.0.0_dp) THEN
          CALL timings_report_print(iw,r_timings,sort_by_self_time,cost_type_energy,para_env)
       ENDIF
       CALL timings_report_print(iw,r_timings,sort_by_self_time,cost_type_time,para_env)
    
       !Write the callgraph, if desired by user
       CALL section_vals_val_get(root_section,"GLOBAL%CALLGRAPH",i_val=cg_mode,error=error)
       IF(cg_mode /= CALLGRAPH_NONE) THEN
          CALL section_vals_val_get(root_section,"GLOBAL%CALLGRAPH_FILE_NAME",c_val=cg_filename,error=error)
          IF(LEN_TRIM(cg_filename) == 0) cg_filename=TRIM(logger%iter_info%project_name)
          IF(cg_mode==CALLGRAPH_ALL)& !incorporate mpi-rank into filename 
             cg_filename = TRIM(cg_filename)//"_"//TRIM(ADJUSTL(cp_to_string(para_env%mepos)))
          IF(iw>0) THEN
             WRITE (UNIT=iw,FMT="(T2,3X,A)") "Writing callgraph to: "//TRIM(cg_filename)//".callgraph"
             WRITE (UNIT=iw,FMT="()") 
             WRITE (UNIT=iw,FMT="(T2,A)") "-------------------------------------------------------------------------------"
          ENDIF
          IF(cg_mode==CALLGRAPH_ALL .OR. para_env%mepos==para_env%source)&
             CALL timings_report_callgraph(TRIM(cg_filename)//".callgraph")
       END IF
       
       CALL cp_print_key_finished_output(iw,logger,root_section,&
            "GLOBAL%TIMINGS",error=error)
    END IF

    
    CALL rm_mp_perf_env()
    CALL rm_timer_env()

    IF (para_env%ionode) THEN
       iw=cp_print_key_unit_nr(logger,root_section,"GLOBAL%PROGRAM_RUN_INFO",&
            extension=".log",error=error)

       ! Deleting (if existing) the external EXIT files
       delete_it = .FALSE.
       INQUIRE (FILE="EXIT",EXIST=delete_it)
       IF(delete_it) THEN
          CALL open_file(file_name="EXIT",unit_number=unit_exit)
          CALL close_file(unit_number=unit_exit,file_status="DELETE")
       END IF

       delete_it = .FALSE.
       INQUIRE (FILE=TRIM(logger%iter_info%project_name)//".EXIT",EXIST=delete_it)
       IF(delete_it) THEN
          CALL open_file(file_name=TRIM(logger%iter_info%project_name)//".EXIT",unit_number=unit_exit)
          CALL close_file(unit_number=unit_exit,file_status="DELETE")
       END IF

       ! update the runtime enviroment variables
       CALL get_runtime_info()

       ! Just a choice, do not print the CP2K footer if there is a failure
       IF (.NOT.failure) THEN
          CALL cp2k_footer(iw,wdir)
          IF (iw>0) CALL m_flush_internal(iw)
       END IF

       CALL cp_print_key_finished_output(iw,logger,root_section,&
            "GLOBAL%PROGRAM_RUN_INFO", error=error)
    END IF
    ! Release message passing environment
    CALL cp_rm_default_logger()

  END SUBROUTINE cp2k_finalize

END MODULE environment
