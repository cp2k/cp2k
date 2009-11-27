!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
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
       compile_arch, compile_date, compile_host, compile_lastcvs, cp2k_home, &
       cp2k_version, cp2k_year, get_runtime_info, id_cp2k_version, &
       r_host_name, r_pid, r_user_name
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
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
       debug_print_level, do_cp2k, do_eip, do_farming, do_fft_acml, &
       do_fft_cu, do_fft_essl, do_fft_fftw, do_fft_fftw2, do_fft_fftw3, &
       do_fft_mkl, do_fft_sci, do_fft_sg, do_fist, do_kg, do_qs, do_test, &
       energy_run, high_print_level, id_development_version, low_print_level, &
       medium_print_level, mol_dyn_run, none_run, silent_print_level
  USE input_cp2k,                      ONLY: create_global_section
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
                                             m_flush_internal,&
                                             m_memory_details
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
  USE reference_manager,               ONLY: print_all_references,&
                                             print_format_journal
  USE string_utilities,                ONLY: ascii_to_string,&
                                             integer_to_string,&
                                             string_to_ascii
  USE termination,                     ONLY: set_error_unit,&
                                             stop_memory,&
                                             stop_program
  USE timings,                         ONLY: add_timer_env,&
                                             rm_timer_env,&
                                             timeprint,&
                                             timeset,&
                                             timestop,&
                                             trace_debug,&
                                             use_HPM
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'environment'

  ! *** Public subroutines ***
  PUBLIC :: cp2k_finalize, cp2k_init, cp2k_read, cp2k_setup

CONTAINS

! *****************************************************************************
!> \brief Initializes a CP2K run (setting of the global environment variables)
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

    CHARACTER(LEN=default_string_length)     :: dev_flag, string
    INTEGER                                  :: ipe, istat, my_output_unit, &
                                                strend, strstart
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: all_pid
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: all_host
    TYPE(cp_logger_type), POINTER            :: logger
#if defined(__HPM)
#include "f_hpm.h"
#endif

    ! create a timer_env
    CALL add_timer_env()

    ! message passing performance
    CALL add_mp_perf_env()

    !   *** Set the default logical output and error unit number ***
    !   *** write to screen in case of error on non-io nodes, and hope this is appropriate
    IF (para_env%ionode) THEN
       CALL set_error_unit(output_unit)
    ELSE
       CALL set_error_unit(6)
    ENDIF

    !   *** Set flag if this is a development version
    dev_flag = ""
    IF (id_cp2k_version==id_development_version) dev_flag=" (Development Version)"

    !   *** Init the default logger
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
    IF (my_output_unit>0) THEN
       DO ipe=1,para_env%num_pe
          CALL ascii_to_string(all_host(:,ipe),string)
          WRITE (UNIT=my_output_unit,FMT="(T2,A,T63,I8,T71,I10)")&
               TRIM(r_user_name)//"@"//TRIM(string)//&
               " has created process number",ipe-1,all_pid(ipe)
       END DO
    END IF
    DEALLOCATE (all_pid,STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,__LINE__,"all_pid")
    END IF
    DEALLOCATE (all_host,STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,__LINE__,"all_host")
    END IF

    !   *** Initialize timing ***
    CALL timeset ("CP2K",globenv%handle)

    !   *** Print header ***
    CALL cp2k_header(my_output_unit,wdir)

    IF (my_output_unit>0) THEN
       WRITE (UNIT=my_output_unit,FMT="(/,T2,A,T31,A50)")&
            "CP2K| version string: ",&
            ADJUSTR(TRIM(cp2k_version)//TRIM(dev_flag))
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
       strend=MAX(1,LEN(compile_lastcvs)-2)
       strstart=MAX(2,strend-50)
       WRITE (UNIT=my_output_unit,FMT="(T2,A,T31,A50)")&
            "CP2K| Last CVS entry",&
            ADJUSTR(compile_lastcvs(strstart:strend))
       WRITE (UNIT=my_output_unit,FMT="(T2,A,T31,A50)")&
            "CP2K| Input file name",&
            ADJUSTR(TRIM(input_file_name))
       CALL m_flush_internal(my_output_unit)
    END IF

  END SUBROUTINE cp2k_init

! *****************************************************************************
!> \brief read part of cp2k_init
!> \param globenv the globenv
!> \note
!>      The following routines need to be synchronized wrt. adding/removing
!>      of the default environments (logging, perormance,error):
!>      environment:cp2k_init, environment:cp2k_finalize,
!>      f77_interface:f_env_add_defaults, f77_interface:f_env_rm_defaults,
!>      f77_interface:create_force_env, f77_interface:destroy_force_env
!> \author fawzi
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
!> \param globenv the global environment to initialize
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      if possible do the initializations here as the environement
!>      (error,...) is setup, instaed of cp2k_init
!> \author fawzi
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
!> \note
!>      Should not be required anymore once everything is converted
!>      to get information directly from the input structure
!> \par History
!>      06-2005 [created]
!> \author MI
! *****************************************************************************
  SUBROUTINE read_global_section(root_section,para_env,globenv,error)

    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_global_section', &
      routineP = moduleN//':'//routineN, start_section_label = "GLOBAL"

    CHARACTER(len=6)                         :: print_level_string
    CHARACTER(len=default_path_length) :: basis_set_file_name, &
      coord_file_name, geminal_file_name, mm_potential_file_name, &
      potential_file_name
    CHARACTER(len=default_string_length)     :: env_num, project_name
    INTEGER :: i_fft, iforce_eval, method_name_id, n_rep_val, nforce_eval, &
      num_threads, output_unit, print_level
    INTEGER(kind=int_8) :: Buffers, Buffers_avr, Buffers_max, Buffers_min, &
      Cached, Cached_avr, Cached_max, Cached_min, MemFree, MemFree_avr, &
      MemFree_max, MemFree_min, MemLikelyFree, MemLikelyFree_avr, &
      MemLikelyFree_max, MemLikelyFree_min, MemTotal, MemTotal_avr, &
      MemTotal_max, MemTotal_min, Slab, Slab_avr, Slab_max, Slab_min, &
      SReclaimable, SReclaimable_avr, SReclaimable_max, SReclaimable_min
    INTEGER, DIMENSION(:), POINTER           :: i_force_eval
    LOGICAL                                  :: ata, efl, explicit, failure, &
                                                trace, trace_master, &
                                                trace_sync
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(enumeration_type), POINTER          :: enum1, enum2
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section
    TYPE(section_vals_type), POINTER         :: dft_section, &
                                                force_env_sections, &
                                                global_section, qmmm_section, &
                                                subsys_section

!$  INTEGER :: omp_get_num_threads

    failure = .FALSE.
    NULLIFY(dft_section,global_section, i_force_eval)

    global_section => section_vals_get_subs_vals(root_section,"GLOBAL",error=error)
    CALL section_vals_val_get(global_section,"BLACS_GRID",i_val=globenv%blacs_grid_layout,error=error)
    CALL section_vals_val_get(global_section,"BLACS_REPEATABLE",l_val=globenv%blacs_repeatable,error=error)
    CALL section_vals_val_get(global_section,"PREFERRED_FFT_LIBRARY",i_val=i_fft,error=error)

    CALL section_vals_val_get(global_section,"PRINT_LEVEL",i_val=print_level,error=error)
    CALL section_vals_val_get(global_section,"PROGRAM_NAME",i_val=globenv%prog_name_id,error=error)
    CALL section_vals_val_get(global_section,"FFT_POOL_SCRATCH_LIMIT",i_val=globenv%fft_pool_scratch_limit,error=error)
    CALL section_vals_val_get(global_section,"FFTW_PLAN_TYPE",i_val=globenv%fftw_plan_type,error=error)
    CALL section_vals_val_get(global_section,"FFTW_ARRAYS_ALIGNED",l_val=globenv%fftw_arrays_aligned,error=error)
    CALL section_vals_val_get(global_section,"PROJECT_NAME",c_val=project_name,error=error)
    CALL section_vals_val_get(global_section,"FFTW_WISDOM_FILE_NAME",c_val=globenv%fftw_wisdom_file_name,error=error)
    CALL section_vals_val_get(global_section,"RUN_TYPE",i_val=globenv%run_type_id,error=error)
    CALL section_vals_val_get(global_section,"WALLTIME",r_val=globenv%cp2k_target_time,error=error)
    CALL section_vals_val_get(global_section,"TRACE",l_val=trace,error=error)
    CALL section_vals_val_get(global_section,"TRACE_SYNC",l_val=trace_sync,error=error)
    CALL section_vals_val_get(global_section,"TRACE_MASTER",l_val=trace_MASTER,error=error)
    CALL section_vals_val_get(global_section,"FLUSH_SHOULD_FLUSH",l_val=flush_should_flush,error=error)
    force_env_sections => section_vals_get_subs_vals(root_section,"FORCE_EVAL",error=error)
    CALL section_vals_get(force_env_sections, n_repetition=nforce_eval, error=error)
    logger=>cp_error_get_logger(error)
    output_unit=cp_print_key_unit_nr(logger,global_section,"PROGRAM_RUN_INFO",&
         extension=".log",error=error)
    IF(trace) CALL trace_debug("start",trace_sync,trace_master,para_env)

    SELECT CASE(i_fft)
    CASE(do_fft_sg)
       globenv%default_fft_library="FFTSG"
    CASE(do_fft_fftw2)
       globenv%default_fft_library="FFTW2"
       CALL cite_reference(Frigo2005)
    CASE(do_fft_fftw3)
       globenv%default_fft_library="FFTW3"
       CALL cite_reference(Frigo2005)
    CASE(do_fft_essl)
       globenv%default_fft_library="FFTESSL"
    CASE(do_fft_acml)
       globenv%default_fft_library="FFTACML"
    CASE(do_fft_mkl)
       globenv%default_fft_library="FFTMKL"
    CASE(do_fft_sci)
       globenv%default_fft_library="FFTSCI"
    CASE(do_fft_cu)
       globenv%default_fft_library="FFTCU"
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,&
            " unknown FFT type",para_env)
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
             CASE (do_kg)
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
            " FARMING program supports only NONE as run type",para_env)
    ENDIF

    IF( globenv%prog_name_id == do_test .AND. globenv%run_type_id /= none_run) &
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


    num_threads=1
!$omp parallel
!$ num_threads = omp_get_num_threads()
!$omp end parallel
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
         "unknown print_level",para_env)
       END SELECT

       WRITE (UNIT=output_unit,FMT="(T2,A,T75,A6)")&
            start_section_label//"| Global print level",print_level_string
       WRITE (UNIT=output_unit,FMT="(T2,A,T75,I6)")&
            start_section_label//"| Total number of message passing processes",&
            para_env%num_pe,&
            start_section_label//"| Number of threads for this process",&
            num_threads,&
            start_section_label//"| This output is from process",para_env%mepos
       WRITE (UNIT=output_unit,FMT='()')

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

    END IF
    CALL cp_print_key_finished_output(output_unit,logger,global_section,&
         "PROGRAM_RUN_INFO", error=error)

  END SUBROUTINE read_global_section

! *****************************************************************************
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

    CHARACTER(LEN=10)                        :: fftw_name
    COMPLEX(KIND=dp), DIMENSION(4, 4, 4)     :: zz
    INTEGER                                  :: output_unit, stat
    INTEGER, DIMENSION(3)                    :: n
    LOGICAL                                  :: try_fftw
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: global_section

    global_section => section_vals_get_subs_vals(root_section,"GLOBAL",error=error)
    CALL read_global_section(root_section,para_env,globenv,error=error)
    logger => cp_error_get_logger(error)
    output_unit=cp_print_key_unit_nr(logger,global_section,"PROGRAM_RUN_INFO",&
         extension=".log",error=error)

    n(:) = 4
    zz(:,:,:) = 0.0_dp

    !
    ! Setup the FFT library
    ! If the user has specified PREFERRED_FFT_LIBRARY try that first (default FFTW)
    ! If that one is not available, try FFTW (unless it has been tried already)
    ! If FFTW is not available use FFTSG
    !
    IF (globenv%default_fft_library.EQ."FFTW2" .OR. globenv%default_fft_library.EQ."FFTW3") THEN
       try_fftw=.FALSE.
    ELSE
       try_fftw=.TRUE.
       SELECT CASE(do_fft_fftw)
       CASE(do_fft_fftw2)
          fftw_name="FFTW2"
       CASE(do_fft_fftw3)
          fftw_name="FFTW3"
       CASE DEFAULT
          CALL stop_program(routineN,moduleN,__LINE__,&
               " unknown FFT type",para_env)
       END SELECT
    ENDIF

    !   *** Initialize FFT library with the user's prefered FFT library ***
    CALL init_fft(fftlib=TRIM(globenv%default_fft_library),&
         alltoall=section_get_lval(global_section,"ALLTOALL_SGL",error), &
         fftsg_sizes=.NOT.section_get_lval(global_section,"EXTENDED_FFT_LENGTHS",error), &
         pool_limit=globenv%fft_pool_scratch_limit,&
         wisdom_file=globenv%fftw_wisdom_file_name,&
         plan_style=globenv%fftw_plan_type,&
         arrays_aligned=globenv%fftw_arrays_aligned,&
         error=error)

    !   *** Check for FFT library ***
    CALL fft3d(1,n,zz,status=stat)
    IF (stat /= 0) THEN

       IF (try_fftw) THEN
          IF (output_unit > 0) THEN
             WRITE(output_unit,'(A,A,T55,A)') &
                  " WARNING : FFT library "//TRIM(globenv%default_fft_library)//&
                  " is not available "," Trying "//TRIM(fftw_name)
          ENDIF
          globenv%default_fft_library=fftw_name
          CALL init_fft(fftlib=TRIM(globenv%default_fft_library),&
               alltoall=section_get_lval(global_section,"ALLTOALL_SGL",error), &
               fftsg_sizes=.NOT.section_get_lval(global_section,"EXTENDED_FFT_LENGTHS",error), &
               pool_limit=globenv%fft_pool_scratch_limit,&
               wisdom_file=globenv%fftw_wisdom_file_name,&
               plan_style=globenv%fftw_plan_type,&
               arrays_aligned=globenv%fftw_arrays_aligned,&
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
               arrays_aligned=globenv%fftw_arrays_aligned,&
               error=error)

          CALL fft3d(1,n,zz,status=stat)
          IF (stat /= 0) THEN
             CALL stop_program(routineN,moduleN,__LINE__,"FFTSG not functional....",para_env)
          ENDIF

       ENDIF

    END IF

    CALL cp_print_key_finished_output(output_unit,logger,global_section,&
         "PROGRAM_RUN_INFO", error=error)

  END SUBROUTINE read_cp2k_section

! *****************************************************************************
!> \brief Writes final timings and banner for CP2K
!> \note
!>      The following routines need to be synchronized wrt. adding/removing
!>      of the default environments (logging, perormance,error):
!>      environment:cp2k_init, environment:cp2k_finalize,
!>      f77_interface:f_env_add_defaults, f77_interface:f_env_rm_defaults,
!>      f77_interface:create_force_env, f77_interface:destroy_force_env
!> \par History
!>      none
!> \author JGH,MK
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
    INTEGER                                  :: iw, unit_exit
    LOGICAL                                  :: delete_it, failure
    REAL(KIND=dp)                            :: r_timings
    TYPE(cp_logger_type), POINTER            :: logger

#if defined(__HPM)
#include "f_hpm.h"
#endif
! look if we inherited a failure, more care is needed if so
! i.e. the input is most likely not available

    failure=.FALSE.
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

    ! finalize the fft (i.e. writes the wisdom)
    CALL finalize_fft(para_env,error)

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

    IF (.NOT. failure) THEN
       iw = cp_print_key_unit_nr(logger,root_section,"GLOBAL%TIMINGS",&
            extension=".Log",error=error)

       r_timings = section_get_rval(root_section,"GLOBAL%TIMINGS%THRESHOLD",error)
       CALL timeprint(iw,r_timings,para_env)

       CALL cp_print_key_finished_output(iw,logger,root_section,&
            "GLOBAL%TIMINGS",error=error)
    ENDIF

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
