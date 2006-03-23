!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2004  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/environment [1.0] *
!!
!!   NAME
!!     environment
!!
!!   FUNCTION
!!     Sets up and terminates the global environment variables
!!
!!   AUTHOR
!!     JGH,MK
!!
!!   MODIFICATION HISTORY
!!     - init_mathcon is now called (27.02.2001,MK)
!!     - Merged with Quickstep MODULE start_program_run (17.01.2002,MK)
!!     - Compile information added (16.01.2002,MK)
!!     - Merged with MODULE cp2k_input, some rearrangements (30.10.2002,MK)
!!
!!******
!******************************************************************************

MODULE environment
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE fft_tools,                       ONLY: fft3d,&
                                             init_fft
  USE gamma,                           ONLY: deallocate_md_ftable
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: &
       do_eip, do_farming, do_fist, do_kg, do_qs, do_tbmd, do_test, &
       energy_run, mol_dyn_run, none_run, program_names, run_type_names
  USE input_section_types,             ONLY: section_get_lval,&
                                             section_get_rval,&
                                             section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp,&
                                             int_size,&
                                             print_kind_info
  USE machine,                         ONLY: m_datum,&
                                             m_getcwd,&
                                             m_getlog,&
                                             m_getpid,&
                                             m_hostnm
  USE message_passing,                 ONLY: add_mp_perf_env,&
                                             describe_mp_perf_env,&
                                             mp_bcast,&
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
  USE string_utilities,                ONLY: ascii_to_string,&
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

#if defined(__COMPILE_ARCH)
  CHARACTER(LEN=*), PARAMETER :: compile_arch = __COMPILE_ARCH
#else
  CHARACTER(LEN=*), PARAMETER :: compile_arch = "unknown: -D__COMPILE_ARCH=?"
#endif

#if defined(__COMPILE_DATE)
  CHARACTER(LEN=*), PARAMETER :: compile_date = __COMPILE_DATE
#else
  CHARACTER(LEN=*), PARAMETER :: compile_date = "unknown: -D__COMPILE_DATE=?"
#endif

#if defined(__COMPILE_HOST)
  CHARACTER(LEN=*), PARAMETER :: compile_host = __COMPILE_HOST
#else
  CHARACTER(LEN=*), PARAMETER :: compile_host = "unknown: -D__COMPILE_HOST=?"
#endif

#if defined(__COMPILE_LASTCVS)
  CHARACTER(LEN=*), PARAMETER :: compile_lastcvs = __COMPILE_LASTCVS
#else
  CHARACTER(LEN=*), PARAMETER :: compile_lastcvs = "/No CVS/Entry in src dir//"
#endif

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'environment'
  INTEGER :: handle = 0

! *** Public subroutines ***

  PUBLIC :: cp2k_finalize,&
            cp2k_init, cp2k_read, cp2k_setup

  PUBLIC :: compile_arch,compile_date,compile_host,compile_lastcvs

!******************************************************************************

CONTAINS

!******************************************************************************
!!****f* environment/cp2k_init [1.0] *
!!
!!   NAME
!!     cp2k_init
!!
!!   FUNCTION
!!     Initializes a CP2K run (setting of the global environment variables)
!!
!!   AUTHOR
!!     JGH,MK
!!
!!   MODIFICATION HISTORY
!!     JGH (28.11.2001) : default for pp_library_path
!!     - print keys added (17.01.2002, MK)
!!     - merged with cp2k_input (30.10.2002,MK)
!!
!!******
!******************************************************************************

  SUBROUTINE cp2k_init(globenv)

    TYPE(global_environment_type), POINTER :: globenv

    CHARACTER(LEN=*), PARAMETER :: routineN = "cp2k_init"

!   *** Local variables ***

    CHARACTER(LEN=default_path_length) :: cwd
    CHARACTER(LEN=30)                  :: host_name,string,user_name
    CHARACTER(LEN=26)                  :: datx
    INTEGER                            :: ipe,istat,l,output_unit,pid,&
                                          strstart,strend
    LOGICAL           :: ionode

    INTEGER, DIMENSION(:), ALLOCATABLE   :: all_pid
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: all_host
    TYPE(cp_logger_type), POINTER :: logger

#if defined(__HPM)
#include "f_hpm.h"    
#endif

!   ---------------------------------------------------------------------------

    CALL m_getpid(pid)
    CALL m_getlog(user_name)
    CALL m_hostnm(host_name)

    ! create a timer_env
    CALL add_timer_env()

    ! message passing performance
    CALL add_mp_perf_env()

!   *** Set the default logical output and error unit number ***
!   *** write to screen in case of error on non-io nodes, and hope this is appropriate
    IF (globenv%ionode) THEN
       CALL set_error_unit(globenv%scr)
    ELSE
       CALL set_error_unit(6)
    ENDIF

!   *** Init the default logger 

    IF (globenv%para_env%source==globenv%para_env%mepos) THEN
       output_unit=globenv%scr
    ELSE
       output_unit=-1
    END IF
    NULLIFY(logger)
    CALL cp_logger_create(logger,para_env=globenv%para_env,&
         default_global_unit_nr=output_unit, &
         close_global_unit_on_dealloc=.FALSE.)
    CALL cp_add_default_logger(logger)
    CALL cp_logger_release(logger)

    ionode = globenv%ionode

    IF (ionode) THEN

      CALL m_datum(datx)
      CALL m_getcwd(cwd)

    END IF

    CALL mp_bcast(globenv%input_file_name,globenv%source,globenv%group)


!   *** Print a list of all started processes ***

    ALLOCATE (all_pid(globenv%num_pe),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineN,moduleN,__LINE__,"all_pid",&
                       globenv%num_pe*int_size)
    END IF
    all_pid(:) = 0
    all_pid(globenv%mepos+1) = pid
    CALL mp_sum(all_pid,globenv%group)
    ALLOCATE (all_host(30,globenv%num_pe),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineN,moduleN,__LINE__,"all_host",&
                       30*globenv%num_pe*int_size)
    END IF
    all_host(:,:) = 0
    CALL string_to_ascii(host_name,all_host(:,globenv%mepos+1))
    CALL mp_sum(all_host,globenv%group)
    IF (globenv%ionode) THEN
      DO ipe=1,globenv%num_pe
        CALL ascii_to_string(all_host(:,ipe),string)
        WRITE (UNIT=globenv%scr,FMT="(T2,A,T63,I8,T71,I10)")&
          TRIM(user_name)//"@"//TRIM(string)//&
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

    CALL timeset ("CP2K","I","",handle)
    output_unit=globenv%scr

    IF (ionode) THEN
      WRITE (UNIT=output_unit,&
             FMT="(/,A,T55,A26,/,A,T51,A30,/,A,T51,A30,/,A,T71,I10)")&
        "  **** **** ******  **  PROGRAM STARTED AT ",ADJUSTR(datx),&
        " ***** ** ***  *** **   PROGRAM STARTED ON ",ADJUSTR(host_name),&
        " **    ****   ******    PROGRAM STARTED BY ",ADJUSTR(user_name),&
        " ***** **    ** ** **   PROGRAM PROCESS ID ",pid
      WRITE (UNIT=output_unit,FMT="(A,T44,A37)")&
        "  **** **  *******  **  PROGRAM STARTED IN ",ADJUSTR(cwd(1:37))
      DO l=38,LEN_TRIM(cwd),37
        WRITE (UNIT=output_unit,FMT="(T44,A)")&
          cwd(l:MIN(LEN_TRIM(cwd),l+36))
      END DO
      WRITE (UNIT=output_unit,FMT="(T2,A,T31,A50)")&
        "CP2K| Program compiled at",&
        ADJUSTR(compile_date(1:MIN(50,LEN(compile_date))))
      WRITE (UNIT=output_unit,FMT="(T2,A,T31,A50)")&
        "CP2K| Program compiled on",&
        ADJUSTR(compile_host(1:MIN(50,LEN(compile_host))))
      WRITE (UNIT=output_unit,FMT="(T2,A,T31,A50)")&
        "CP2K| Program compiled for",&
        ADJUSTR(compile_arch(1:MIN(50,LEN(compile_arch))))
      strend=MAX(1,LEN(compile_lastcvs)-2)
      strstart=MAX(2,strend-50)
      WRITE (UNIT=output_unit,FMT="(T2,A,T31,A50)")&
        "CP2K| Last CVS entry",&
        ADJUSTR(compile_lastcvs(strstart:strend))
    END IF


  END SUBROUTINE cp2k_init

! *****************************************************************************

!!****f* environment/cp2k_read [1.0] *
!!
!!   NAME
!!     cp2k_read
!!
!!   FUNCTION
!!     read part of cp2k_init
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - globenv: the globenv
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp2k_read(globenv,error)
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp2k_read', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: c_val
    INTEGER                                  :: iw
    TYPE(cp_logger_type), POINTER            :: logger

!   *** Read the input/output section ***

    logger => cp_get_default_logger()
#if defined(__HPM)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         globenv%input_file,"GLOBAL%PERFORMANCE/HPM",error=error),cp_p_file)) THEN
        CALL f_hpminit(globenv%para_env%mepos,"CP2K") 
        use_HPM = .TRUE.
    ENDIF
#endif

!   *** Read the CP2K section ***

    CALL read_cp2k_section(globenv)

    ! try to use better names for the local log if it is not too late
    CALL section_vals_val_get(globenv%input_file,"GLOBAL%OUTPUT_FILE_NAME",&
         c_val=c_val,error=error)
    IF (c_val/="") THEN
       CALL cp_logger_set(logger,&
            local_filename=TRIM(c_val)//"_localLog")
    END IF
    CALL section_vals_val_get(globenv%input_file,"GLOBAL%PROJECT",&
         c_val=c_val,error=error)
    IF (c_val/="") THEN
       CALL cp_logger_set(logger,&
            local_filename=TRIM(c_val)//"_localLog")
    END IF
    logger%iter_info%project_name=c_val
    logger%iter_info%print_level=globenv%print_level

    iw = cp_print_key_unit_nr(logger,globenv%input_file,"GLOBAL%PRINT/BASIC_DATA_TYPES",&
         extension=".Log",error=error)
    IF (iw>0) CALL print_kind_info(iw)
    CALL cp_print_key_finished_output(iw,logger,globenv%input_file,&
         "GLOBAL%PRINT/BASIC_DATA_TYPES",error=error)

    iw = cp_print_key_unit_nr(logger,globenv%input_file,"GLOBAL%PRINT/PHYSCON",&
         extension=".Log",error=error)
    IF (iw>0) CALL write_physcon(iw)
    CALL cp_print_key_finished_output(iw,logger,globenv%input_file,&
         "GLOBAL%PRINT/PHYSCON",error=error)    

  END SUBROUTINE cp2k_read
!***************************************************************************

!!****f* environment/cp2k_setup [1.0] *
!!
!!   NAME
!!     cp2k_setup
!!
!!   FUNCTION
!!     globenv initializations that need the input and error
!!
!!   NOTES
!!     if possible do the initializations here as the environement
!!     (error,...) is setup, instaed of cp2k_init
!!
!!   INPUTS
!!     - globenv: the global environment to initialize
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************

  SUBROUTINE cp2k_setup(globenv,error)

    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp2k_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: input_seed, iw
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(3, 2)           :: initial_seed
    TYPE(cp_logger_type), POINTER            :: logger

! -------------------------------------------------------------------------

    failure = .FALSE.

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    IF (.NOT.failure) THEN

      ! Initialize the parallel random number generator

      CALL init_rng()
      iw = cp_print_key_unit_nr(logger,globenv%input_file,"GLOBAL%PRINT/RNG_MATRICES",&
                                extension=".Log",error=error)
      IF (iw > 0) THEN
        CALL write_rng_matrices(iw)
      END IF
      CALL cp_print_key_finished_output(iw,logger,globenv%input_file,&
                                        "GLOBAL%PRINT/RNG_MATRICES",&
                                        error=error)

      ! Initialize a global normally Gaussian distributed (pseudo)random number stream

      CALL section_vals_val_get(globenv%input_file,"GLOBAL%SEED",i_val=input_seed)
      initial_seed(:,:) = REAL(input_seed,KIND=dp)
      CALL create_rng_stream(rng_stream=globenv%gaussian_rng_stream,&
                             name="Global Gaussian random numbers",&
                             distribution_type=GAUSSIAN,&
                             seed=initial_seed,&
                             extended_precision=.TRUE.,&
                             error=error)

      iw = cp_print_key_unit_nr(logger,globenv%input_file,"GLOBAL%PRINT/RNG_CHECK",&
                                extension=".Log",error=error)
      IF (iw > 0) THEN
        CALL check_rng(iw,globenv%ionode,error)
      END IF
      CALL cp_print_key_finished_output(iw,logger,globenv%input_file,&
                                        "GLOBAL%PRINT/RNG_CHECK",&
                                        error=error)

      iw = cp_print_key_unit_nr(logger,globenv%input_file,"GLOBAL%PRINT/GLOBAL_GAUSSIAN_RNG",&
                                extension=".Log",error=error)
      IF (iw > 0) THEN
        CALL write_rng_stream(globenv%gaussian_rng_stream,iw,write_all=.TRUE.,error=error)
      END IF
      CALL cp_print_key_finished_output(iw,logger,globenv%input_file,&
                                        "GLOBAL%PRINT/GLOBAL_GAUSSIAN_RNG",&
                                        error=error)

    END IF

  END SUBROUTINE cp2k_setup

!***************************************************************************

!!****f* environment/read_global_section [1.0] *
!!
!!   NAME
!!     read_global_section
!!
!!   FUNCTION
!!     read the global sectionof new input
!!
!!   NOTES
!!     Should not be required anymore once everything is converted 
!!     to get information directly from the input structure
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     06-2005 [created]
!!
!!******
!******************************************************************************
  SUBROUTINE read_global_section(globenv,error)

    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_global_section', &
      start_section_label = "GLOBAL"

    CHARACTER(len=default_string_length) :: basis_set_file_name, &
      mm_potential_file_name, potential_file_name, restart_file_name
    INTEGER                                  :: num_threads, output_unit
    LOGICAL                                  :: ata, failure, trace
    TYPE(section_vals_type), POINTER         :: dft_section, global_section

!$  INTEGER :: omp_get_num_threads

    failure = .FALSE.
    trace = .FALSE.
    NULLIFY(dft_section,global_section)

    global_section => section_vals_get_subs_vals(globenv%input_file,"GLOBAL",error=error)

    CALL section_vals_val_get(global_section,"BLACS_GRID",i_val=globenv%blacs_grid_layout,error=error)
    CALL section_vals_val_get(global_section,"FFT_LIBRARY",c_val=globenv%default_fft_library,error=error)
    CALL section_vals_val_get(global_section,"PRINT_LEVEL",i_val=globenv%print_level,error=error)
    CALL section_vals_val_get(global_section,"PROGRAM_NAME",i_val=globenv%prog_name_id,error=error)
    CALL section_vals_val_get(global_section,"PROJECT_NAME",c_val=globenv%project_name,error=error)
    CALL section_vals_val_get(global_section,"RUN_TYPE",i_val=globenv%run_type_id,error=error)
    CALL section_vals_val_get(global_section,"WALLTIME",r_val=globenv%cp2k_target_time,error=error)
    CALL section_vals_val_get(global_section,"TRACE",l_val=trace,error=error)
    IF(trace) CALL trace_debug("start",groupid=globenv%group)

    globenv%coord_file_name = globenv%input_file_name

    dft_section => section_vals_get_subs_vals(globenv%input_file,"FORCE_EVAL%DFT",error=error)

    CALL section_vals_val_get(dft_section,"BASIS_SET_FILE_NAME",&
         c_val=basis_set_file_name , ignore_required=.TRUE., error=error)
    CALL section_vals_val_get(dft_section,"POTENTIAL_FILE_NAME",&
         c_val=potential_file_name , ignore_required=.TRUE., error=error)
    CALL section_vals_val_get(dft_section,"RESTART_FILE_NAME",&
         c_val=restart_file_name , ignore_required=.TRUE., error=error)
    CALL section_vals_val_get(globenv%input_file,"FORCE_EVAL%QMMM%MM_POTENTIAL_FILE_NAME",&
         c_val=mm_potential_file_name, ignore_required=.TRUE., error=error)

    IF (globenv%run_type_id == 0) THEN
       SELECT CASE (globenv%prog_name_id)
       CASE (do_farming)
          globenv%run_type_id = none_run
       CASE (do_fist)
          globenv%run_type_id = mol_dyn_run
       CASE (do_kg)
          globenv%run_type_id = mol_dyn_run
       CASE (do_tbmd)
          globenv%run_type_id = mol_dyn_run
       CASE (do_eip)
          globenv%run_type_id = mol_dyn_run
       CASE (do_qs)
          globenv%run_type_id = energy_run
       CASE (do_test)
          globenv%run_type_id = none_run
       END SELECT
    END IF

    IF(globenv%prog_name_id == do_farming .AND. globenv%run_type_id /= none_run) &
       CALL stop_program(routineN,moduleN,__LINE__,&
               " FARMING program supports only NONE as run type",globenv)

    IF(globenv%prog_name_id == do_test .AND. globenv%run_type_id /= none_run) &
       CALL stop_program(routineN,moduleN,__LINE__,&
               "TEST program supports only NONE as run type",globenv)

    num_threads=1
!$omp parallel
!$ num_threads = omp_get_num_threads()
!$omp end parallel
    IF (globenv%ionode.AND.(globenv%print_level > 0)) THEN
      output_unit = globenv%scr
      WRITE (UNIT=output_unit,FMT=*)
      WRITE (UNIT=output_unit,FMT="(T2,A,T41,A)")&
        start_section_label//"| Input file name",&
        ADJUSTR(globenv%input_file_name(:40)),&
        start_section_label//"| Basis set file name",&
        ADJUSTR(basis_set_file_name(:40)),&
        start_section_label//"| Potential file name",&
        ADJUSTR(potential_file_name(:40)),&
        start_section_label//"| MM Potential file name",&
        ADJUSTR(mm_potential_file_name(:40)),&
        start_section_label//"| Restart file name",&
        ADJUSTR(restart_file_name(:40)),&
        start_section_label//"| Coordinate file name",&
        ADJUSTR(globenv%coord_file_name(:40)),&
        start_section_label//"| Program name",&
        ADJUSTR(program_names(globenv%prog_name_id)(:40)),&
        start_section_label//"| Project name",&
        ADJUSTR(globenv%project_name(:40)),&
        start_section_label//"| Default FFT library",&
        ADJUSTR(globenv%default_fft_library(:40)),&
        start_section_label//"| Run type",&
        ADJUSTR(run_type_names(globenv%run_type_id)(:40))

      CALL section_vals_val_get(global_section,"ALLTOALL_SGL",l_val=ata,error=error)
      WRITE (UNIT=output_unit,FMT="(T2,A,T80,L1)")&
        start_section_label//"| All-to-all communication in single precision",ata

      WRITE (UNIT=output_unit,FMT="(T2,A,T75,I6)")&
        start_section_label//"| Global print level",globenv%print_level,&
        start_section_label//"| Total number of message passing processes",&
                             globenv%num_pe,&
        start_section_label//"| Number of threads for this process",&
                             num_threads,&
        start_section_label//"| This output is from process",globenv%mepos
      WRITE (UNIT=output_unit,FMT='()')
    END IF

  END SUBROUTINE read_global_section

!******************************************************************************
!!****f* environment/read_cp2k_section [1.0] *
!!
!!   NAME
!!     read_cp2k_section
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH,MK
!!
!!   MODIFICATION HISTORY
!!     2-Dec-2000 (JGH) added default fft library
!!
!!******
!******************************************************************************

  SUBROUTINE read_cp2k_section(globenv,error)

    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'read_cp2k_section'&
                                                , start_section_label = "CP2K"

    COMPLEX(KIND=dp), DIMENSION(4, 4, 4)     :: zz
    INTEGER                                  :: stat
    INTEGER, DIMENSION(3)                    :: n
    LOGICAL                                  :: was_present
    TYPE(section_vals_type), POINTER         :: global_section

!   ---------------------------------------------------------------------------

    was_present = .FALSE.

    global_section => section_vals_get_subs_vals(globenv%input_file,"GLOBAL",error=error)
    CALL section_vals_get(global_section,explicit=was_present,error=error)
    CALL read_global_section(globenv)

!   *** Initialize FFT library ***
    CALL init_fft(fftlib=TRIM(globenv%default_fft_library),&
                  alltoall=section_get_lval(global_section,"ALLTOALL_SGL",error))

!   *** Check for FFT library ***
    n(:) = 4
    zz(:,:,:) = 0.0_dp
    CALL fft3d(1,n,zz,status=stat)
    IF (stat /= 0) THEN
      IF (globenv%ionode) THEN
         WRITE(globenv%scr,'(A,A,T55,A)') &
           " WARNING : FFT library "//TRIM(globenv%default_fft_library)//&
           " is not available "," Trying FFTSG as a default"
      ENDIF
      globenv%default_fft_library="FFTSG"
      CALL init_fft(fftlib=TRIM(globenv%default_fft_library),&
                    alltoall=section_get_lval(global_section,"ALLTOALL_SGL",error))
      CALL fft3d(1,n,zz,status=stat) 
      IF (stat /= 0) THEN
           CALL stop_program(routineN,moduleN,__LINE__,&
                             "The default FFT library <"//&
                             TRIM(globenv%default_fft_library)//&
                             "> is not available",globenv)
      ENDIF
    END IF

  END SUBROUTINE read_cp2k_section

!******************************************************************************
!!****f* environment/cp2k_finalize [1.0] *
!!
!!   NAME
!!     cp2k_finalize
!!
!!   FUNCTION
!!     Writes final timings and banner for CP2K
!!
!!   AUTHOR
!!     JGH,MK
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!******
!******************************************************************************

  SUBROUTINE cp2k_finalize(globenv,error)

    TYPE(global_environment_type), POINTER :: globenv
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    CHARACTER(LEN=default_path_length) :: cwd
    CHARACTER(LEN=30)                  :: host_name,user_name
    CHARACTER(LEN=26)                  :: datx
    INTEGER                            :: l,iw,pid,unit_exit
    LOGICAL                            :: delete_it,failure
    TYPE(cp_logger_type), POINTER      :: logger
    REAL(KIND=dp)                      :: r_timings
#if defined(__HPM)
#include "f_hpm.h"
#endif

    ! -------------------------------------------------------------------------

    ! look if we inherited a failure, more care is needed if so 
    ! i.e. the input is most likely not available

    failure=.FALSE.
    CALL cp_error_check(error,failure) 

    ! Clean up

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    CALL deallocate_spherical_harmonics()
    CALL deallocate_orbital_pointers()
    CALL deallocate_md_ftable()

    ! Write message passing performance info

    IF (.NOT. failure) THEN
       CALL describe_mp_perf_env ( globenv%group,globenv%scr )
    ENDIF
    CALL rm_mp_perf_env()

#if defined(__HPM)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         globenv%input_file,"GLOBAL%PERFORMANCE/HPM",error=error),cp_p_file)) THEN
       CALL f_hpmterminate(globenv%group) 
       use_HPM=.FALSE.
    ENDIF
#endif

    CALL timestop(0.0_dp,handle) ! corresponding the "CP2K" in cp2k_init

    IF (.NOT. failure) THEN
       iw = cp_print_key_unit_nr(logger,globenv%input_file,"GLOBAL%TIMINGS",&
            extension=".Log",error=error)
       IF (iw>0) THEN
          r_timings = section_get_rval(globenv%input_file,"GLOBAL%TIMINGS%THRESHOLD",error)
          CALL timeprint(iw,r_timings)
       END IF
       CALL cp_print_key_finished_output(iw,logger,globenv%input_file,&
            "GLOBAL%TIMINGS",error=error)
    ENDIF

    CALL rm_timer_env()

    IF (globenv%ionode) THEN

      iw = globenv%scr

      CALL m_datum(datx)
      CALL m_getcwd(cwd)
      CALL m_getlog(user_name)
      CALL m_getpid(pid)
      CALL m_hostnm(host_name)

      delete_it = .FALSE.
      INQUIRE (FILE="EXIT",EXIST=delete_it)
      IF(delete_it) THEN
         CALL open_file(file_name="EXIT",unit_number=unit_exit)
         CALL close_file(unit_number=unit_exit,file_status="DELETE")
      END IF

      ! Just a choice, do not print the final banner if there is a failure

      IF (globenv%ionode .AND. .NOT. failure) THEN
        WRITE (UNIT=iw,&
               FMT="(A,T55,A26,/,A,T51,A30,/,A,T51,A30,/,A,T71,I10)")&
          "  **** **** ******  **  PROGRAM ENDED AT   ",ADJUSTR(datx),&
          " ***** ** ***  *** **   PROGRAM RAN ON     ",ADJUSTR(host_name),&
          " **    ****   ******    PROGRAM RAN BY     ",ADJUSTR(user_name),&
          " ***** **    ** ** **   PROGRAM PROCESS ID ",pid
        WRITE (UNIT=iw,FMT="(A,T44,A37)")&
          "  **** **  *******  **  PROGRAM STARTED IN ",ADJUSTR(cwd(1:37))
        DO l=38,LEN_TRIM(cwd),37
          WRITE (UNIT=iw,FMT="(T44,A)") cwd(l:MIN(LEN_TRIM(cwd),l+36))
        END DO
      END IF

    END IF

    ! Release message passing environment

    CALL cp_rm_default_logger()

  END SUBROUTINE cp2k_finalize

  ! ***************************************************************************

END MODULE environment
