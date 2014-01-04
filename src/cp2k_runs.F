!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
MODULE cp2k_runs

  USE atom,                            ONLY: atom_code
  USE bsse,                            ONLY: do_bsse_calculation
  USE cell_opt,                        ONLY: cp_cell_opt
  USE cp2k_debug,                      ONLY: cp2k_debug_energy_and_forces
  USE cp2k_info,                       ONLY: compile_date,&
                                             compile_revision,&
                                             cp2k_version,&
                                             cp2k_year,&
                                             enable_color_tags
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_config,&
                                             cp_dbcsr_finalize_lib,&
                                             cp_dbcsr_init_lib,&
                                             cp_dbcsr_print_config
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_unit_nr,&
                                             cp_rm_iter_level
  USE cp_para_env,                     ONLY: cp_para_env_create,&
                                             cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_search_string
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE cp_units,                        ONLY: cp_unit_set_create,&
                                             cp_unit_set_release,&
                                             cp_unit_set_type,&
                                             print_all_units
  USE cuda_device,                     ONLY: cuda_device_init
  USE cuda_memory,                     ONLY: cuda_device_mem_init,&
                                             cuda_device_mem_release
  USE environment,                     ONLY: cp2k_finalize,&
                                             cp2k_init,&
                                             cp2k_read,&
                                             cp2k_setup
  USE f77_blas
  USE f77_interface,                   ONLY: create_force_env,&
                                             destroy_force_env,&
                                             f77_default_para_env => default_para_env,&
                                             f_env_add_defaults,&
                                             f_env_rm_defaults,&
                                             f_env_type
  USE farming_methods,                 ONLY: do_deadlock,&
                                             do_nothing,&
                                             do_wait,&
                                             farming_parse_input,&
                                             get_next_job
  USE farming_types,                   ONLY: deallocate_farming_env,&
                                             farming_env_type,&
                                             init_farming_env,&
                                             job_finished,&
                                             job_running
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE geo_opt,                         ONLY: cp_geo_opt
  USE global_types,                    ONLY: global_environment_type,&
                                             globenv_create,&
                                             globenv_release,&
                                             globenv_retain
  USE input_constants,                 ONLY: &
       bsse_run, cell_opt_run, debug_run, do_atom, do_band, do_cp2k, do_ep, &
       do_farming, do_fist, do_mixed, do_opt_basis, do_optimize_input, &
       do_qmmm, do_qs, do_swarm, do_tamc, do_test, do_tree_mc, &
       do_tree_mc_ana, driver_run, ehrenfest, electronic_spectra_run, energy_force_run, &
       energy_run, geo_opt_run, linear_response_run, mol_dyn_run, &
       mon_car_run, none_run, pint_run, real_time_propagation, tree_mc_run, &
       vib_anal
  USE input_cp2k,                      ONLY: create_cp2k_input_reading,&
                                             create_cp2k_root_section,&
                                             create_global_section,&
                                             empty_initial_variables
  USE input_cp2k_check,                ONLY: check_cp2k_input
  USE input_keyword_types,             ONLY: keyword_release
  USE input_parsing,                   ONLY: section_vals_parse
  USE input_section_types,             ONLY: &
       section_describe_html, section_describe_index_html, section_release, &
       section_type, section_typo_match, section_vals_create, &
       section_vals_get_subs_vals, section_vals_release, section_vals_retain, &
       section_vals_type, section_vals_val_get, section_vals_write, &
       write_section_xml
  USE ipi_driver,                      ONLY: run_driver
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE library_tests,                   ONLY: lib_test
  USE machine,                         ONLY: m_chdir,&
                                             m_flush,&
                                             m_getcwd,&
                                             m_walltime
  USE mc_run,                          ONLY: do_mon_car
  USE md_run,                          ONLY: qs_mol_dyn
  USE message_passing,                 ONLY: &
       mp_any_source, mp_bcast, mp_comm_dup, mp_comm_free, mp_comm_split, &
       mp_environ, mp_recv, mp_send, mp_sum, mp_sync
  USE neb_methods,                     ONLY: neb
  USE optimize_basis,                  ONLY: run_optimize_basis
  USE optimize_input,                  ONLY: run_optimize_input
  USE pint_methods,                    ONLY: do_pint_run
  USE qs_linres_module,                ONLY: linres_calculation
  USE qs_tddfpt_module,                ONLY: tddfpt_calculation
  USE reference_manager,               ONLY: print_all_references,&
                                             print_format_html
  USE rt_propagation,                  ONLY: rt_prop_setup
  USE string_utilities,                ONLY: html_entity_table,&
                                             write_html_tables
  USE swarm,                           ONLY: run_swarm
  USE tamc_run,                        ONLY: qs_tamc
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE tmc_setup,                       ONLY: do_analyze_files,&
                                             do_tmc
  USE vibrational_analysis,            ONLY: vb_anal
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC  :: cp2k_run, write_cp2k_html_manual, write_xml_file, run_input, &
             farming_run
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp2k_runs'
CONTAINS

! *****************************************************************************
!> \brief performs an instance of a cp2k run
!>      cp2k_run_OK == .TRUE. if everything went fine
!> \param input_file_name name of the file to be opened for input
!> \param ouput_unit unit to which output should be written
!> \param mpi_env the mpi environement for the run (the caller is
!>                        responsible to get rid of it)
!> \note
!>      para_env should be a valid communicator
!>      output_unit should be writeable by at least the lowest rank of the mpi group
!>
!>      recursive because a given run_type might need to be able to perform
!>      another cp2k_run as part of its job (e.g. farming, classical equilibration, ...)
!>
!>      the idea is that a cp2k instance should be able to run with just three
!>      arguments, i.e. a given input file, output unit, mpi communicator.
!>      giving these three to cp2k_run should produce a valid run.
!>      the only task of the PROGRAM cp2k is to create valid instances of the
!>      above arguments. Ideally, anything that is called afterwards should be
!>      able to run simultaneously / multithreaded / sequential / parallel / ...
!>      and able to fail safe
!> \author Joost VandeVondele
! *****************************************************************************
  RECURSIVE FUNCTION cp2k_run(input_file_name,output_unit,mpi_comm) RESULT(cp2k_run_OK)

    CHARACTER(LEN=*), INTENT(IN)             :: input_file_name
    INTEGER, INTENT(IN)                      :: output_unit, mpi_comm
    LOGICAL                                  :: cp2k_run_OK

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp2k_run', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=512)                       :: matching_string(15), &
                                                unknown_string
    INTEGER :: f_env_handle, i, ierr, iter_level, matching_rank(15), &
      method_name_id, new_env_id, prog_name_id, run_type_id, unit_nr
    LOGICAL                                  :: echo_input, failure, &
                                                I_was_ionode, was_present
    TYPE(cp_error_type)                      :: error, suberror
    TYPE(cp_logger_type), POINTER            :: logger, sublogger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_type), POINTER              :: input_structure
    TYPE(section_vals_type), POINTER         :: glob_section, input_file, &
                                                root_section

    cp2k_run_OK = .FALSE.
    failure=.FALSE.
    was_present = .FALSE.
    NULLIFY(para_env, f_env)
    CALL cp_para_env_create(para_env, group=mpi_comm,&
         owns_group=.FALSE.,error=error)

    CALL cp_dbcsr_init_lib (group=mpi_comm, error=error)

    CALL cuda_device_init (group=mpi_comm, error=error)

    NULLIFY(globenv, force_env)

    CALL cp_error_init(error, stop_level=cp_fatal_level)

    ! parse the input
    input_file => create_cp2k_input_reading(input_file_name,initial_variables=empty_initial_variables, &
                                            para_env=para_env,error=error)

    ! put failure to true if this is what happened on parsing the input
    CALL cp_error_check(error,failure)
    IF (failure) THEN
       logger => cp_error_get_logger(error)
       unit_nr=cp_logger_get_default_io_unit(logger)
       IF (unit_nr>0) THEN
          CALL cp_error_get(error,info=unknown_string)
          WRITE(unit_nr,'(T2,A)') ""
          WRITE(unit_nr,'(T2,A)') "Looking for words in the input similar to the unknown: "
          WRITE(unit_nr,'(T4,A)') ''''//TRIM(unknown_string)//''''
          WRITE(unit_nr,'(T2,A)') ""
          CALL m_flush(unit_nr)
          matching_rank=0
          matching_string=""
          NULLIFY(input_structure)
          CALL cp_error_init(suberror,template_error=error)
          CALL create_cp2k_root_section(input_structure,error=suberror)
          CALL cp_error_dealloc_ref(suberror)
          CALL section_typo_match(input_structure,unknown_string,location_string="", &
                                  matching_rank=matching_rank,matching_string=matching_string,error=error)
          CALL section_release(input_structure,error=error)
          DO I=1,SIZE(matching_rank,1)
             IF (matching_rank(i)>0) THEN
                 WRITE(unit_nr,'(T2,A)') TRIM(matching_string(I))
             ENDIF
          ENDDO
          WRITE(unit_nr,'(T15,A)') ""
          WRITE(unit_nr,'(T15,A)') "CP2K failed to parse the input file."
          WRITE(unit_nr,'(T15,A)') "A full description of the input for this CP2K version"
          WRITE(unit_nr,'(T15,A)') "can be generated using:"
          WRITE(unit_nr,'(T15,A)') ""
          WRITE(unit_nr,'(T15,A)') "cp2k.sopt --html-manual"
          WRITE(unit_nr,'(T15,A)') ""
          WRITE(unit_nr,'(T15,A)') "The manual for the latest version of CP2K is online available:"
          WRITE(unit_nr,'(T15,A)') ""
          WRITE(unit_nr,'(T15,A)') "http://manual.cp2k.org/trunk"
          WRITE(unit_nr,'(T15,A)') ""
          WRITE(unit_nr,'(T15,A)') "If this input was an input of a previous version"
          WRITE(unit_nr,'(T15,A)') "of CP2K, you can try to convert it with --permissive-echo."
          WRITE(unit_nr,'(T15,A)') "However, this will just ignore the unknown keywords ..."
          WRITE(unit_nr,'(T15,A)') ""
          CALL m_flush(unit_nr)
       ENDIF
       CALL mp_sync(para_env%group)
    ELSE
       glob_section => section_vals_get_subs_vals(input_file,"GLOBAL",error=error)
       CALL section_vals_val_get(glob_section,"ECHO_INPUT",l_val=echo_input,&
            error=error)
       logger => cp_error_get_logger(error)
       IF (echo_input) THEN
          CALL section_vals_write(input_file,&
               unit_nr=cp_logger_get_default_io_unit(logger),&
               hide_root=.TRUE., hide_defaults=.FALSE., error=error)
       END IF

       IF (.NOT.failure) THEN
          CALL check_cp2k_input(input_file,para_env=para_env,output_unit=output_unit,error=error)
       END IF
       CALL cp_error_check(error,failure)

       IF (failure) THEN
          logger => cp_error_get_logger(error)
          unit_nr=cp_logger_get_default_io_unit(logger)
          IF (unit_nr>0) THEN
             WRITE(unit_nr,'(T15,A)') ""
             WRITE(unit_nr,'(T15,A)') "CP2K failed to parse the input."
             WRITE(unit_nr,'(T15,A)') "It was presumably generated by an incompatible CP2K version."
             WRITE(unit_nr,'(T15,A)') "try --permissive-echo to generate a compatible input file (to be checked!!!)"
             WRITE(unit_nr,'(T15,A)') ""
             CALL m_flush(unit_nr)
          ENDIF
       ELSE
          root_section=>input_file
          CALL section_vals_val_get(input_file,"GLOBAL%PROGRAM_NAME",&
               i_val=prog_name_id,error=error)
          CALL section_vals_val_get(input_file,"GLOBAL%RUN_TYPE",&
               i_val=run_type_id,error=error)

          IF (prog_name_id/=do_cp2k) THEN
             ! initial setup (cp2k does in in the creation of the force_env)
             CALL globenv_create(globenv, error=error)
             ! XXXXXXXXX
             ! root_section => input_file
             ! XXXXXXXXX
             CALL section_vals_retain(input_file,error=error)
             CALL cp2k_init(para_env, output_unit, globenv, input_file_name=input_file_name)
             CALL cp_error_init(suberror)
             CALL cp2k_read(root_section,para_env,globenv,error=error)
             CALL cp2k_setup(root_section,para_env, globenv,suberror)
          END IF

          CALL cp_dbcsr_config (root_section, error)
          IF (output_unit > 0 .AND.&
               cp_logger_would_log(logger, cp_note_level)) THEN
             CALL cp_dbcsr_print_config (unit_nr=output_unit, error=error)
             WRITE (UNIT=output_unit,FMT='()')
          ENDIF

          IF (prog_name_id/=do_farming) THEN
             ! if PW_CUDA is enabled we need to get cudaStreams and
             ! cudaEvents on the GPU
             CALL cuda_device_mem_init(root_section, error)
          ENDIF


          SELECT CASE (prog_name_id)
          CASE (do_atom)
             globenv%run_type_id = none_run
             CALL atom_code(root_section,suberror)
          CASE (do_optimize_input)
             CALL run_optimize_input(root_section, para_env, suberror)
          CASE (do_swarm)
             CALL run_swarm(root_section, para_env, globenv, input_file_name, suberror)
          CASE (do_farming)
             ! Hack: DBCSR should be uninitialized when entering farming.
             ! But, we don't want to change the public f77_interface.
             ! TODO: refactor cp2k's startup code
             CALL cp_dbcsr_finalize_lib(group=mpi_comm, output_unit=0, error=error)
             CALL farming_run ( root_section, para_env, suberror )
             CALL cp_dbcsr_init_lib(group=mpi_comm, error=error)
          CASE (do_opt_basis)
             CALL run_optimize_basis ( root_section, para_env, suberror )
             globenv%run_type_id = none_run
          CASE(do_cp2k)
             CALL create_force_env(new_env_id,input_path=input_file_name,&
                  output_path="__STD_OUT__",mpi_comm=para_env%group,&
                  output_unit=output_unit,&
                  owns_out_unit=.FALSE.,&
                  input=input_file,ierr=ierr)
             CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
             CALL f_env_add_defaults(new_env_id,f_env,new_error=suberror,&
                  failure=failure,handle=f_env_handle)
             force_env => f_env%force_env
             CALL force_env_get(force_env,globenv=globenv,error=suberror)
             CALL globenv_retain(globenv,error=suberror)
             CALL section_vals_val_get(force_env%force_env_section,"METHOD",i_val=method_name_id,error=error)

          CASE (do_test)
             CALL lib_test(root_section,para_env,globenv,suberror)
          CASE (do_tree_mc) ! TMC entry point
             CALL do_tmc ( root_section, para_env, globenv, error=suberror )
          CASE (do_tree_mc_ana)
             CALL do_analyze_files(root_section, para_env, globenv, error=suberror )
          CASE default
             CPAssert(.FALSE.,cp_failure_level,routineP,suberror,failure)
          END SELECT
          CALL section_vals_release(input_file,error=suberror)
          SELECT CASE (globenv%run_type_id)
          CASE (pint_run)
             CALL do_pint_run( para_env, root_section, globenv, error=suberror )
          CASE (none_run, tree_mc_run)
             ! do nothing
          CASE (energy_run, energy_force_run)
             IF(  method_name_id /= do_qs .AND.&
                  method_name_id /= do_qmmm .AND.&
                  method_name_id /= do_ep .AND.&
                  method_name_id /= do_mixed .AND.&
                  method_name_id /= do_fist )  &
                  CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                  "Energy/Force run not available for all methods "//&
                  CPSourceFileRef,&
                  suberror,failure)

             sublogger => cp_error_get_logger(error)
             CALL cp_add_iter_level(sublogger%iter_info,"JUST_ENERGY",&
                  n_rlevel_new=iter_level,&
                  error=suberror)
             SELECT CASE(globenv%run_type_id)
             CASE (energy_run)
                CALL force_env_calc_energy_force(force_env,calc_force=.FALSE.,error=suberror)
             CASE(energy_force_run)
                CALL force_env_calc_energy_force(force_env,calc_force=.TRUE. ,error=suberror)
             CASE default
                CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
             END SELECT
             CALL cp_rm_iter_level(sublogger%iter_info,level_name="JUST_ENERGY",n_rlevel_att=iter_level,&
                  error=suberror)
          CASE (mol_dyn_run)
             CALL qs_mol_dyn ( force_env, globenv, error=suberror )
          CASE (geo_opt_run)
             CALL cp_geo_opt(force_env,globenv,error=suberror)
          CASE (cell_opt_run)
             CALL cp_cell_opt(force_env,globenv,error=suberror)
          CASE (mon_car_run)
             CALL do_mon_car ( force_env, globenv, input_file_name, error=suberror )
          CASE (do_tamc)
             CALL qs_tamc ( force_env, globenv, input_file_name, error=suberror )
          CASE (electronic_spectra_run)
             IF(method_name_id /= do_qs) &
                  CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                  "Electron spectra available only with Quickstep. "//&
                  CPSourceFileRef,&
                  suberror,failure)
             CALL force_env_calc_energy_force(force_env,calc_force=.FALSE.,&
                  error=suberror)
             CALL tddfpt_calculation(force_env%qs_env, error=error)
          CASE (real_time_propagation)
             IF(method_name_id /= do_qs) &
                  CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                  "Real time propagation needs METHOD QS. "//&
                  CPSourceFileRef,&
                  suberror,failure)
             force_env%qs_env%dft_control%rtp_control%fixed_ions=.TRUE.
             CALL rt_prop_setup(force_env,error)
          CASE (ehrenfest)
             IF(method_name_id /= do_qs) &
                  CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                  "Ehrenfest dynamics needs METHOD QS "//&
                  CPSourceFileRef,&
                  suberror,failure)
             force_env%qs_env%dft_control%rtp_control%fixed_ions=.FALSE.
             CALL qs_mol_dyn ( force_env, globenv, error=suberror )
          CASE (bsse_run)
             CALL do_bsse_calculation(force_env, globenv, error=suberror)
          CASE (linear_response_run)
             IF(method_name_id /= do_qs .AND. &
                method_name_id /= do_qmmm) &
                CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                     "Property calculations by Linear Response only within the QS or QMMM program "//&
                     CPSourceFileRef,&
                     suberror,failure)
             ! The Ground State is needed, it can be read from Restart
             CALL force_env_calc_energy_force(force_env,calc_force=.FALSE.,linres=.TRUE.,&
                                              error=suberror)
             CALL linres_calculation(force_env, error=suberror)
          CASE (debug_run)
             SELECT CASE(method_name_id)
             CASE(do_qs, do_qmmm, do_ep, do_fist)
                CALL cp2k_debug_energy_and_forces(force_env, error=suberror)
             CASE DEFAULT
                CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                     "Debug run available only with QS, FIST, and QMMM program "//&
                     CPSourceFileRef,&
                     suberror,failure)
             END SELECT
          CASE (vib_anal)
             CALL vb_anal(root_section,para_env,globenv,error=suberror)
          CASE (do_band)
             CALL neb(root_section,para_env,globenv,error=suberror)
          CASE default
             CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT


          ! if CUDA is enabled we need to free memory on the graphics card
          CALL cuda_device_mem_release()
          CALL cp_dbcsr_finalize_lib(group=mpi_comm, output_unit=output_unit, error=error)


          CALL cp_error_propagate_error(suberror,routineP,&
               "cp2k ending run", error=error,failure=failure)
          IF (prog_name_id==do_cp2k) THEN
             f_env%force_env => force_env ! for mc
             CALL globenv_retain(globenv,error=error)!mc
             CALL globenv_release(force_env%globenv,error=error) !mc
             force_env%globenv => globenv !mc
             CALL f_env_rm_defaults(f_env,error=suberror,ierr=ierr,&
                  handle=f_env_handle)
             CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
             CALL destroy_force_env(new_env_id,ierr=ierr)
             CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
             CALL globenv_release(globenv,error=error)
          ELSE
             CALL cp_error_dealloc_ref(suberror)
             I_was_ionode=para_env%ionode
             CALL cp2k_finalize(root_section,para_env,globenv,error=error)
             CPPostconditionNoFail(globenv%ref_count==1,cp_failure_level,routineP,error)
             CALL section_vals_release(root_section,error=error)
             CALL globenv_release(globenv,error=error)
          END IF

          ! if PW_CUDA is enabled we need to free cudaStreams, cudaEvents
          ! and stored FFT-plans on the GPU
          CALL cuda_device_mem_release()

       END IF

    END IF

    CALL cp_error_dealloc_ref(error)
    CALL cp_para_env_release(para_env,error=error)
    cp2k_run_OK = .NOT. failure

  END FUNCTION cp2k_run

! *****************************************************************************
!> \brief performs a farming run that performs several independent cp2k_runs
!> \note
!>      needs to be part of this module as the cp2k_run -> farming_run -> cp2k_run
!>      calling style creates a hard circular dependency
!> \author Joost VandeVondele
! *****************************************************************************
  RECURSIVE SUBROUTINE farming_run(root_section,para_env, error)
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'farming_run', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: slave_status_done = -3, &
                                                slave_status_wait = -4

    CHARACTER(len=7)                         :: label
    CHARACTER(LEN=default_path_length)       :: output_file
    CHARACTER(LEN=default_string_length)     :: str
    INTEGER :: dest, handle, i, i_job_to_restart, ierr, ijob, ijob_current, &
      ijob_end, ijob_start, iunit, n_jobs_to_run, new_group, new_output_unit, &
      new_rank, new_size, ngroups, num_slaves, output_unit, primus_slave, &
      slave_group, slave_rank, source, stat, tag, todo
    INTEGER, DIMENSION(:), POINTER           :: group_distribution, &
                                                master_slave_partition, &
                                                slave_distribution, &
                                                slave_status
    LOGICAL                                  :: failure, found, master, &
                                                run_OK, slave
    REAL(KIND=dp)                            :: t1, t2
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: waittime
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_parser_type), POINTER            :: my_parser
    TYPE(cp_unit_set_type), POINTER          :: default_units
    TYPE(farming_env_type), POINTER          :: farming_env
    TYPE(section_type), POINTER              :: g_section
    TYPE(section_vals_type), POINTER         :: g_data


    failure=.FALSE.
    ! the primus of all slaves, talks to the master on topics concerning all slaves
    CALL timeset(routineN,handle)
    NULLIFY(my_parser,g_section,g_data,default_units)

    logger => cp_error_get_logger(error)
    output_unit=cp_print_key_unit_nr(logger,root_section,"FARMING%PROGRAM_RUN_INFO",&
                                         extension=".log",error=error)

    IF (output_unit>0) WRITE(output_unit,FMT="(T2,A)") "FARMING| Hi, welcome on this farm!"

    ALLOCATE(farming_env,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL init_farming_env(farming_env)
    ! remember where we started
    CALL m_getcwd(farming_env%cwd)
    CALL farming_parse_input(farming_env,root_section,para_env,error)

    ! the full mpi group is first split in a slave group and a master group, the latter being at most 1 process
    slave=.TRUE.
    master=.FALSE.
    IF (farming_env%master_slave) THEN
       IF (output_unit>0) WRITE(output_unit,FMT="(T2,A)") "FARMING| using a master-slave setup"

       ALLOCATE(master_slave_partition(0:1))
       master_slave_partition=(/1,para_env%num_pe-1/)
       ALLOCATE(group_distribution(0:para_env%num_pe-1),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       CALL mp_comm_split(para_env%group,slave_group,ngroups,group_distribution,&
                   n_subgroups=2, group_partition=master_slave_partition)
       DEALLOCATE(master_slave_partition)
       DEALLOCATE(group_distribution)
       CALL mp_environ(num_slaves,slave_rank,slave_group)

       IF (para_env%mepos==0) THEN
           slave =.FALSE.
           master=.TRUE.
           ! on the master node, num_slaves corresponds to the size of the master group
           ! due to the mp_environ call.
           CPPostcondition(num_slaves==1,cp_failure_level,routineP,error,failure)
           num_slaves=para_env%num_pe-1
           slave_rank=-1
       ENDIF
       CPPostcondition(num_slaves==para_env%num_pe-1,cp_failure_level,routineP,error,failure)
    ELSE
       ! all processes are slaves
       IF (output_unit>0) WRITE(output_unit,FMT="(T2,A)") "FARMING| using a slave-only setup"
       CALL mp_comm_dup(para_env%group,slave_group)
       CALL mp_environ(num_slaves,slave_rank,slave_group)
    ENDIF
    IF (output_unit>0) WRITE(output_unit,FMT="(T2,A,I0)") "FARMING| number of slaves ",num_slaves

    ! keep track of which para_env rank is which slave/master
    ALLOCATE(slave_distribution(0:para_env%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    slave_distribution=0
    slave_distribution(para_env%mepos)=slave_rank
    CALL mp_sum(slave_distribution,para_env%group)
    ! we do have a primus inter pares
    primus_slave=0
    DO i=1,para_env%num_pe-1
       IF (slave_distribution(i)==0) primus_slave=i
    ENDDO

    ! split the current communicator for the slaves
    ! in a new_group, new_size and new_rank according to the number of groups required according to the input
    ALLOCATE(group_distribution(0:num_slaves-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    group_distribution=-1
    IF (slave) THEN
        IF (farming_env%group_size_wish_set) THEN
           farming_env%group_size_wish=MIN(farming_env%group_size_wish,para_env%num_pe)
           CALL mp_comm_split(slave_group,new_group,ngroups,group_distribution,&
                   subgroup_min_size=farming_env%group_size_wish)
        ELSE IF (farming_env%ngroup_wish_set) THEN
           CALL mp_comm_split(slave_group,new_group,ngroups,group_distribution,&
                   n_subgroups=farming_env%ngroup_wish,&
                   group_partition=farming_env%group_partition)
        ELSE
           CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
        ENDIF
        CALL mp_environ(new_size,new_rank,new_group)
    ENDIF

    ! transfer the info about the slave group distribution to the master
    IF (farming_env%master_slave) THEN
       IF (para_env%mepos==primus_slave) THEN
           tag=1
           CALL mp_send(group_distribution,0,tag,para_env%group)
           tag=2
           CALL mp_send(ngroups,0,tag,para_env%group)
       ENDIF
       IF (para_env%mepos==0) THEN
           tag=1
           CALL mp_recv(group_distribution,primus_slave,tag,para_env%group)
           tag=2
           CALL mp_recv(ngroups,primus_slave,tag,para_env%group)
       ENDIF
    ENDIF

    ! write info on group distribution
    IF (output_unit>0) THEN
        WRITE(output_unit,FMT="(T2,A,T71,I10)") "FARMING| Number of created MPI (slave) groups:",ngroups
        WRITE(output_unit,FMT="(T2,A)",ADVANCE="NO") "FARMING| MPI (slave) process to group correspondence:"
        DO i=0,num_slaves-1
           IF (MODULO(i,4)==0) WRITE(output_unit,*)
           WRITE(output_unit,FMT='(A3,I6,A3,I6,A1)',ADVANCE="NO")&
                "  (",i," : ",group_distribution(i),")"
        END DO
        WRITE(output_unit,*)
    ENDIF

    ! protect about too many jobs being run in single go. Not more jobs are allowed than the number in the input file
    ! and determine the future restart point
    IF (farming_env%cycle) THEN
       n_jobs_to_run=farming_env%max_steps*ngroups
       i_job_to_restart=MODULO(farming_env%restart_n+n_jobs_to_run-1,farming_env%njobs)+1
    ELSE
       n_jobs_to_run=MIN(farming_env%njobs,farming_env%max_steps*ngroups)
       n_jobs_to_run=MIN(n_jobs_to_run,farming_env%njobs-farming_env%restart_n+1)
       i_job_to_restart=n_jobs_to_run+farming_env%restart_n
    ENDIF

    ! and write the restart now, that's the point where the next job starts, even if this one is running
    iunit=cp_print_key_unit_nr(logger,root_section,"FARMING%RESTART",&
         extension=".restart",error=error)
    IF (iunit>0) THEN
       WRITE(iunit,*) i_job_to_restart
    ENDIF
    CALL cp_print_key_finished_output(iunit,logger,root_section,"FARMING%RESTART",error=error)


    ! this is the job range to be executed.
    ijob_start=farming_env%restart_n
    ijob_end=ijob_start+n_jobs_to_run-1
    IF (output_unit>0 .AND. ijob_end-ijob_start<0) THEN
       WRITE(output_unit,FMT="(T2,A)") "FARMING| --- WARNING --- NO JOBS NEED EXECUTION ? "
       WRITE(output_unit,FMT="(T2,A)") "FARMING| is the cycle keyword required ?"
       WRITE(output_unit,FMT="(T2,A)") "FARMING| or is a stray RESTART file present ?"
       WRITE(output_unit,FMT="(T2,A)") "FARMING| or is the group_size requested smaller than the number of CPUs?"
    ENDIF

    ! actual executions of the jobs in two different modes
    IF (farming_env%master_slave) THEN
       IF (slave) THEN
          ! keep on doing work until master has decided otherwise
          todo=do_wait
          DO
             IF (new_rank==0) THEN
                ! the head slave tells the master he's done or ready to start
                ! the message tells what has been done lately
                tag=1
                dest=0
                CALL mp_send(todo,dest,tag,para_env%group)

                ! gets the new todo item
                tag=2
                source=0
                CALL mp_recv(todo,source,tag,para_env%group)

                ! and informs his peer slaves
                CALL mp_bcast(todo,0,new_group)
             ELSE
                CALL mp_bcast(todo,0,new_group)
             ENDIF

             ! if the todo is do_nothing we are flagged to quit. Otherwise it is the job number
             SELECT CASE(todo)
             CASE(do_wait,do_deadlock)
                ! go for a next round, but we first wait a bit
                t1=m_walltime()
                DO
                  t2=m_walltime()
                  IF (t2-t1>farming_env%wait_time) EXIT
                ENDDO
             CASE(do_nothing)
                EXIT
             CASE(1:)
                CALL execute_job(todo)
             END SELECT
          ENDDO
       ELSE ! master
          ALLOCATE(slave_status(0:ngroups-1),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          slave_status=slave_status_wait
          ijob_current=ijob_start-1

          DO
             IF (ALL(slave_status==slave_status_done)) EXIT

             ! who's the next slave waiting for work
             tag=1
             source=mp_any_source
             CALL mp_recv(todo,source,tag,para_env%group) ! updates source
             IF (todo>0) THEN
                farming_env%Job(todo)%status=job_finished
                IF (output_unit>0) THEN
                   WRITE(output_unit,FMT=*) "Job finished: ",todo
                   CALL m_flush(output_unit)
                ENDIF
             ENDIF

             ! get the next job in line, this could be do_nothing, if we're finished
             CALL get_next_job(farming_env,ijob_start,ijob_end,ijob_current,todo)
             dest=source
             tag =2
             CALL mp_send(todo,dest,tag,para_env%group)

             IF (todo>0) THEN
               farming_env%Job(todo)%status=job_running
               IF (output_unit>0) THEN
                 WRITE(output_unit,FMT=*) "Job: ",todo," Dir: ",TRIM(farming_env%Job(todo)%cwd), &
                                        " assigned to group ",group_distribution(slave_distribution(dest))
                 CALL m_flush(output_unit)
               ENDIF
             ELSE
               IF (todo==do_nothing) THEN
                   slave_status(group_distribution(slave_distribution(dest)))=slave_status_done
                   IF (output_unit>0) THEN
                      WRITE(output_unit,FMT=*) "group done: ",group_distribution(slave_distribution(dest))
                      CALL m_flush(output_unit)
                   ENDIF
               ENDIF
               IF (todo==do_deadlock) THEN
                   IF (output_unit>0) THEN
                      WRITE(output_unit,FMT=*) ""
                      WRITE(output_unit,FMT=*) "FARMING JOB DEADLOCKED ... CIRCULAR DEPENDENCIES"
                      WRITE(output_unit,FMT=*) ""
                      CALL m_flush(output_unit)
                   ENDIF
                   CPPostcondition(todo.NE.do_deadlock,cp_failure_level,routineP,error,failure)
               ENDIF
             ENDIF

          ENDDO

          DEALLOCATE(slave_status,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       ENDIF
    ELSE
       ! this is the non-master-slave mode way of executing the jobs
       ! the i-th job in the input is always executed by the MODULO(i-1,ngroups)-th group
       ! (needed for cyclic runs, we don't want two groups working on the same job)
       IF (output_unit>0) THEN
         IF (ijob_end-ijob_start>=0) THEN
           WRITE(output_unit,FMT="(T2,A)") "FARMING| List of jobs : "
           DO ijob=ijob_start,ijob_end
            i=MODULO(ijob-1,farming_env%njobs)+1
            WRITE(output_unit,FMT=*) "Job: ",i," Dir: ",TRIM(farming_env%Job(i)%cwd)," Input: ", &
              TRIM(farming_env%Job(i)%input)," MPI group:", MODULO(i-1,ngroups)
           ENDDO
         ENDIF
       ENDIF

       DO ijob=ijob_start,ijob_end
          i=MODULO(ijob-1,farming_env%njobs)+1
          ! this farms out the jobs
          IF (MODULO(i-1,ngroups)==group_distribution(slave_rank)) THEN
             IF (output_unit > 0) WRITE(output_unit,FMT="(T2,A,I5.5,A)",ADVANCE="NO") " Running Job ",i, &
                           " in "//TRIM(farming_env%Job(i)%cwd)//"."
             CALL execute_job(i)
             IF (output_unit > 0) THEN
                WRITE(output_unit,FMT="(A)") " Done, output in "//TRIM(output_file)
                CALL m_flush(output_unit)
             ENDIF
          ENDIF
       ENDDO
    ENDIF

    ! keep information about how long each process has to wait
    ! i.e. the load imbalance
    t1=m_walltime()
    CALL mp_sync(para_env%group)
    t2=m_walltime()
    ALLOCATE(waittime(0:para_env%num_pe-1))
    waittime=0.0_dp
    waittime(para_env%mepos)=t2-t1
    CALL mp_sum(waittime,para_env%group)
    IF (output_unit>0) THEN
       WRITE(output_unit,'(T2,A)') "Process idle times [s] at the end of the run"
       DO i=0,para_env%num_pe-1
         WRITE(output_unit,FMT='(A2,I6,A3,F8.3,A1)',ADVANCE="NO")&
                " (",i," : ",waittime(i),")"
         IF (MOD(i+1,4)==0) WRITE(output_unit,'(A)') ""
       ENDDO
       CALL m_flush(output_unit)
    ENDIF
    DEALLOCATE(waittime)

    ! give back the communicators of the split groups
    IF (slave) CALL mp_comm_free(new_group)
    CALL mp_comm_free(slave_group)

    ! and message passing deallocate structures
    DEALLOCATE(group_distribution,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(slave_distribution,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! clean the farming env
    CALL deallocate_farming_env(farming_env)

    CALL cp_print_key_finished_output(output_unit,logger,root_section,&
               "FARMING%PROGRAM_RUN_INFO", error=error)

    CALL timestop(handle)

  CONTAINS
! *****************************************************************************
    SUBROUTINE execute_job(i)
    INTEGER                                  :: i

! change to the new working directory

       CALL m_chdir(TRIM(farming_env%Job(i)%cwd),ierr)
       CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

       ! generate a fresh call to cp2k_run
       IF (new_rank == 0) THEN

          IF (farming_env%Job(i)%output=="") THEN
             ! generate the output file
             WRITE(output_file,'(A12,I5.5)') "FARMING_OUT_",i
             CALL parser_create(my_parser,file_name=TRIM(farming_env%Job(i)%input),error=error)
             label="&GLOBAL"
             CALL parser_search_string(my_parser,label,ignore_case=.TRUE.,found=found,error=error)
             IF (found) THEN
                CALL create_global_section(g_section,error=error)
                CALL section_vals_create(g_data,g_section,error=error)
                CALL cp_unit_set_create(default_units, "OUTPUT",error=error)
                CALL section_vals_parse(g_data,my_parser,default_units,&
                     error=error)
                CALL cp_unit_set_release(default_units,error=error)
                CALL section_vals_val_get(g_data,"PROJECT",&
                     c_val=str, error=error)
                IF (str.NE."") output_file=TRIM(str)//".out"
                CALL section_vals_val_get(g_data,"OUTPUT_FILE_NAME",&
                        c_val=str,error=error)
                IF (str.NE."") output_file=str
                CALL section_vals_release(g_data,error=error)
                CALL section_release(g_section,error=error)
             END IF
             CALL parser_release(my_parser,error=error)
          ELSE
             output_file=farming_env%Job(i)%output
          ENDIF

          CALL open_file(file_name=TRIM(output_file),&
                 file_action="WRITE",&
                 file_status="UNKNOWN",&
                 file_position="APPEND",&
                 unit_number=new_output_unit)
       ELSE
         ! this unit should be negative, otherwise all processors that get a default unit
         ! start writing output (to the same file, adding to confusion).
         ! error handling should be careful, asking for a local output unit if required
         new_output_unit=-1
       ENDIF

       run_OK=cp2k_run(TRIM(farming_env%Job(i)%input),new_output_unit,new_group)

       IF (new_rank == 0) CALL close_file(unit_number=new_output_unit)

       ! change to the original working directory
       CALL m_chdir(TRIM(farming_env%cwd),ierr)
       CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

     END SUBROUTINE execute_job
  END SUBROUTINE farming_run

! *****************************************************************************
!> \brief writes a small html description of the cp2k input
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author joost [fawzi]
! *****************************************************************************
  SUBROUTINE write_cp2k_html_manual(error)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_cp2k_html_manual', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: root_section

    failure=.FALSE.

    IF (.NOT. failure) THEN
       NULLIFY(root_section)
       CALL create_cp2k_root_section(root_section,error=error)
       ! remove the default keyword that ignores things outside the section
       CALL keyword_release(root_section%keywords(0)%keyword,error=error)

       CALL section_describe_html(root_section,"InputReference",0,6,error=error)

       CALL open_file(unit_number=unit_nr,file_name="index.html",&
            file_action="WRITE", file_status="REPLACE")
       WRITE(unit_nr,FMT='(A)') "<HTML><BODY><HEAD><TITLE>The cp2k input structure</TITLE>"
       WRITE(unit_nr,FMT='(A)') "<H1>CP2K input reference</H1>"
       WRITE(unit_nr,FMT='(A)') "<H2>Version information</H2>"
       WRITE(unit_nr,FMT='(A)') &
            "This html manual can be generated automatically from a given cp2k executable "// &
            "using the --html-manual command line option. The manual describes exactly that version of the code. "//&
            "This version has been generated using a cp2k version compiled on "//TRIM(compile_date)//". "//&
            "The SVN source code revision number is "//TRIM(compile_revision)//"."

       WRITE(unit_nr,FMT='(A)') '<H2>Search the online input reference manual</H2>'
       WRITE(unit_nr,FMT='(A)') '<TABLE><TR><TD>'
       ! Yahoo seems to actually search our pages
       WRITE(unit_nr,FMT='(A)') ' <!-- START CODE --> '
       WRITE(unit_nr,FMT='(A)') '         <div id="ysrchForm" style="border:1px solid #7E9DB9;'
       WRITE(unit_nr,FMT='(A)') '                                background:#FFFFFF;'
       WRITE(unit_nr,FMT='(A)') '                                width:300px;'
       WRITE(unit_nr,FMT='(A)') '                                margin:0 auto;'
       WRITE(unit_nr,FMT='(A)') '                                padding:20px; '
       WRITE(unit_nr,FMT='(A)') '                                position:relative;'
       WRITE(unit_nr,FMT='(A)') '                              ">'
       WRITE(unit_nr,FMT='(A)') &
            '           <form id="searchBoxForm_undefined" action="http://search.yahoo.com/search" style="padding:0;">'
       WRITE(unit_nr,FMT='(A)') '                         '
       WRITE(unit_nr,FMT='(A)') '                         <input name="ei" value="UTF-8" type="hidden">'
       WRITE(unit_nr,FMT='(A)') &
            '                         <input name="fr" value="ystg" type="hidden"><div style="padding:0 80px 0 0;zoom:1;">'
       WRITE(unit_nr,FMT='(A)') '                                         <input type="text" id="searchTerm"'
       WRITE(unit_nr,FMT='(A)') '                                                 onFocus="this.style.background=''#fff'';"'
       WRITE(unit_nr,FMT='(A)') &
            '                                                 onBlur="if(this.value=='''')this.style.background='//&
            '''#fff url(http://us.i1.yimg.com/us.yimg.com/i/us/sch/gr/horiz_pwrlogo_red2.gif) 3px center  no-repeat''"'
       WRITE(unit_nr,FMT='(A)') &
            ' name="p" style=" margin:1px 0; width:100%; border:1px solid #7E9DB9; color:#666666; height:18px; '//&
            'padding:0px 3px; background:#fff url(http://us.i1.yimg.com/us.yimg.com/i/us/sch/gr/horiz_pwrlogo_red2.gif)'//&
            ' 3px center no-repeat; position:relative;">'
       WRITE(unit_nr,FMT='(A)') ' <input type="submit" id="btn_undefined" value="Search" '
       WRITE(unit_nr,FMT='(A)') &
            ' style=" padding-bottom:2px; position:absolute; right:20px; top:20px; margin:0px; height:22px; width:65px; ">'
       WRITE(unit_nr,FMT='(A)') ' </div><ul style="color:#666666;'
       WRITE(unit_nr,FMT='(A)') '            font:11px/11px normal Arial, Helvetica, sans-serif;'
       WRITE(unit_nr,FMT='(A)') '            margin:0;'
       WRITE(unit_nr,FMT='(A)') '            padding:0;'
       WRITE(unit_nr,FMT='(A)') '            text-align:left;'
       WRITE(unit_nr,FMT='(A)') '            list-style-type:none;radios"><li style="display:inline;padding-right:10px;">'
       WRITE(unit_nr,FMT='(A)') &
            ' <input name="vs" id="web_radio" value=""    type="radio" style="vertical-align:middle;margin-right:5px; ">'
       WRITE(unit_nr,FMT='(A)') ' <label for="web_radio" style="vertical-align:middle;">Web</label>'
       WRITE(unit_nr,FMT='(A)') ' </li><li style="display:inline;padding-right:10px;">'
       WRITE(unit_nr,FMT='(A)') &
            ' <input name="vs" id="site_radio" value="www.cp2k.org"   checked="checked" type="radio"'//&
            ' style="vertical-align:middle;margin-right:5px; ">'
       WRITE(unit_nr,FMT='(A)') ' <label for="site_radio" style="vertical-align:middle;">this Site</label>'
       WRITE(unit_nr,FMT='(A)') ' </li></ul></form></div>'
       WRITE(unit_nr,FMT='(A)') ' <!-- END CODE --></TD><TD> '

       WRITE(unit_nr,FMT='(A)') '</TD></TR></TABLE>'

       WRITE(unit_nr,FMT='(A)') '<H2>Journal papers</H2>'
       WRITE(unit_nr,FMT='(A)') '<A HREF="references.html">List of references</A> cited in the CP2K input manual.'

       WRITE(unit_nr,FMT='(A)') '<H2>Units of Measurement</H2>'
       WRITE(unit_nr,FMT='(A)') '<A HREF="units.html">List of unites</A> organized into groups based on the physical quantity '//&
            "they measure. The units can be used to specify an alternative unit of measurement for keyword values, for which a "//&
            "default unit has been explicitly defined."

       IF (enable_color_tags) THEN
          WRITE(unit_nr,FMT='(A)') "<H2>Color convention</H2>"
          WRITE(unit_nr,FMT='(A)') &
               "Sections and Keywords supported in the current released version of cp2k are marked "//&
               "with green color. Keywords or Sections specifically designed for the development version or still "//&
               "in a development status (not ready for the released version) are marked with black/blue color."
       END IF

       WRITE(unit_nr,FMT='(A)') "<H2>Internal Input Preprocessor</H2>"
       WRITE(unit_nr,FMT='(A)') "Before the input is parsed, the input is run through a very simple internal preprocessor."
       WRITE(unit_nr,FMT='(A)') "The preprocessor recognizes the following directives independent of capitalization:<BR><DL>"
       WRITE(unit_nr,FMT='(A)') "<DT><B>@INCLUDE 'filename.inc'</B></DT>"
       WRITE(unit_nr,FMT='(A)') "<DD>The file referenced by <I>filename.inc</I> is included into the input file and parsed."
       WRITE(unit_nr,FMT='(A)') "Recursive inclusions are not allowed and the files have to exist in the current working "
       WRITE(unit_nr,FMT='(A)') "directory. There can be only one @INCLUDE statement per line. Single or double quotes "
       WRITE(unit_nr,FMT='(A)') "can be used and <B>have</B> to be used if the filename contains blanks.</DD>"
       WRITE(unit_nr,FMT='(A)') "<DT><B>@SET VAR value</B></DT>"
       WRITE(unit_nr,FMT='(A)') "<DD>Assigns the text <I>value</I> to the preprocessing variable <I>VAR</I>. <I>value</I> "
       WRITE(unit_nr,FMT='(A)') "is the text following <I>VAR</I> with the outer whitespace removed. The variable can be "
       WRITE(unit_nr,FMT='(A)') "recalled with a <I>${VAR}</I> statement. There can be only one @SET statement per line.</DD>"
       WRITE(unit_nr,FMT='(A)') "<DT><B>${VAR}</B></DT>"
       WRITE(unit_nr,FMT='(A)') "<DD>Expand the variable <I>VAR</I>. The text <I>${VAR}</I> is replaced with the value assigned "
       WRITE(unit_nr,FMT='(A)') "to <I>VAR</I> in the last corresponding @SET directive. There can be multiple variables per line. "
       WRITE(unit_nr,FMT='(A)') "The expansion process is repeated until no more variables are found.</DD>"
       WRITE(unit_nr,FMT='(A)') "<DT><B>@IF / @ENDIF</B></DT>"
       WRITE(unit_nr,FMT='(A)') "<DD>Conditional block. The text from the @IF line up to the next line with a valid "
       WRITE(unit_nr,FMT='(A)') "@ENDIF is skipped, if the expression following @IF resolves to <I>false</I>. "
       WRITE(unit_nr,FMT='(A)') "Available expressions are lexical comparisons for equality '==' or inequality '/='."
       WRITE(unit_nr,FMT='(A)') "If none of the two operators are found, a '0' or whitespace resolves to <I>false</I> "
       WRITE(unit_nr,FMT='(A)') "while any text resolves to <I>true</I>. @IF/@ENDIF blocks cannot be nested and "
       WRITE(unit_nr,FMT='(A)') "cannot span across files. There can be only one test (== or /=) per @IF statement.</DD>"
       WRITE(unit_nr,FMT='(A)') "</DL><P>"

       WRITE(unit_nr,FMT='(A)') "<H2>Input structure</H2>"
       WRITE(unit_nr,FMT='(A)') "All sections that can be part of a cp2k input file are shown with their allowed nestings. "
       WRITE(unit_nr,FMT='(A)') "A description of each section, and a list of keywords can be obtained clicking on the links. "
       WRITE(unit_nr,FMT='(A)') "<BR><UL>"

       CALL section_describe_index_html(root_section,"InputReference",unit_nr,error=error)

       WRITE(unit_nr,FMT='(A)') '</UL><BR><hr>Back to the <A HREF="http://www.cp2k.org/">CP2K homepage</A>'
       WRITE(unit_nr,FMT='(A)') "</BODY></HTML>"

       CALL close_file(unit_nr)
       CALL section_release(root_section,error=error)

       ! References
       CALL open_file(unit_number=unit_nr,file_name="references.html",&
            file_action="WRITE", file_status="REPLACE")
       WRITE(unit_nr,FMT='(A)') "<HTML><BODY><HEAD><TITLE>The cp2k literature list</TITLE>"
       WRITE(unit_nr,FMT='(A)') "<H1>CP2K references</H1>"
       CALL print_all_references(sorted=.TRUE.,cited_only=.FALSE.,&
                          FORMAT=print_format_html,unit=unit_nr)
       WRITE(unit_nr,FMT='(A)') '<BR><hr>Back to the <A HREF="http://www.cp2k.org/">CP2K homepage</A>'
       WRITE(unit_nr,FMT='(A)') "</BODY></HTML>"
       CALL close_file(unit_nr)

       ! Units
       CALL open_file(unit_number=unit_nr,file_name="units.html",&
            file_action="WRITE", file_status="REPLACE")
       WRITE(unit_nr,FMT='(A)') "<HTML><BODY><HEAD><TITLE>The cp2k units list</TITLE>"
       WRITE(unit_nr,FMT='(A)') "<H1>CP2K Available Units of Measurement</H1>"
       CALL print_all_units(unit_nr=unit_nr)
       WRITE(unit_nr,FMT='(A)') '<BR><hr>Back to the <A HREF="http://www.cp2k.org/">CP2K homepage</A>'
       WRITE(unit_nr,FMT='(A)') "</BODY></HTML>"
       CALL close_file(unit_nr)
    END IF
  END SUBROUTINE write_cp2k_html_manual

! *****************************************************************************
  SUBROUTINE write_xml_file(error)

    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_xml_file', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ie, is, unit_number
    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: root_section

    failure = .FALSE.
    IF (.NOT.failure) THEN
      NULLIFY(root_section)
      CALL create_cp2k_root_section(root_section,error=error)
      CALL keyword_release(root_section%keywords(0)%keyword,error)
      CALL open_file(unit_number=unit_number,&
                     file_name="cp2k_input.xml",&
                     file_action="WRITE",&
                     file_status="REPLACE")

      WRITE (UNIT=unit_number,FMT="(A)")&
        "<?xml version=""1.0"" encoding=""ISO-8859-1""?>",&
        "<?xml-stylesheet type=""text/xsl"" href=""cp2k_input.xsl""?>"

      !MK Write a HTML translation table
      !MK In principle this is only required for non-standard HTML entities
      is = 0
      ie = 0
      WRITE (UNIT=unit_number,FMT="(A)")&
        "<!DOCTYPE documentElement["
      DO i=1,SIZE(html_entity_table),2
         is = INDEX(html_entity_table(i),"&") + 1
         CPPostcondition((is > 0),cp_failure_level,routineP,error,failure)
         ie = INDEX(html_entity_table(i),";") - 1
         CPPostcondition((ie >= is),cp_failure_level,routineP,error,failure)
         WRITE (UNIT=unit_number,FMT="(A)")&
          "<!ENTITY "//html_entity_table(i)(is:ie)//" """//TRIM(html_entity_table(i+1))//""">"
      END DO
      WRITE (UNIT=unit_number,FMT="(A)")&
        "]>"

      !MK CP2K input structure
      WRITE (UNIT=unit_number,FMT="(A)")&
        "<CP2K_INPUT>",&
        " <CP2K_VERSION>"//TRIM(cp2k_version)//"</CP2K_VERSION>",&
        " <CP2K_YEAR>"//TRIM(cp2k_year)//"</CP2K_YEAR>",&
        " <COMPILE_DATE>"//TRIM(compile_date)//"</COMPILE_DATE>",&
        " <COMPILE_REVISION>"//TRIM(compile_revision)//"</COMPILE_REVISION>"
      DO i=1,root_section%n_subsections
         CALL write_section_xml(root_section%subsections(i)%section,1,unit_number,error)
      END DO
      ! Append HTML entity and tag tables
      CALL write_html_tables(unit_number)
      WRITE (UNIT=unit_number,FMT="(A)") "</CP2K_INPUT>"
      CALL close_file(unit_number=unit_number)
      CALL section_release(root_section,error=error)

      ! References
      CALL open_file(unit_number=unit_number,file_name="references.html",&
                     file_action="WRITE",file_status="REPLACE")
      WRITE(unit_number,FMT='(A)') "<HTML><BODY><HEAD><TITLE>The cp2k literature list</TITLE>"
      WRITE(unit_number,FMT='(A)') "<H1>CP2K references</H1>"
      CALL print_all_references(sorted=.TRUE.,cited_only=.FALSE., &
                                FORMAT=print_format_html,unit=unit_number)
      WRITE(unit_number,FMT='(A)') "</BODY></HTML>"
      CALL close_file(unit_number=unit_number)

      ! Units
      CALL open_file(unit_number=unit_number,file_name="units.html",&
                     file_action="WRITE",file_status="REPLACE")
      WRITE(unit_number,FMT='(A)') "<HTML><BODY><HEAD><TITLE>The cp2k units list</TITLE>"
      WRITE(unit_number,FMT='(A)') "<H1>CP2K Available Units of Measurement</H1>"
      CALL print_all_units(unit_nr=unit_number)
      WRITE(unit_number,FMT='(A)') "</BODY></HTML>"
      CALL close_file(unit_number=unit_number)

    END IF

  END SUBROUTINE write_xml_file

! *****************************************************************************
!> \brief runs the given input
!> \param input_file_path the path of the input file
!> \param output_file_path path of the output file (to which it is appended)
!>        if it is "__STD_OUT__" the unit 6 is used
!> \param mpi_comm the mpi communicator to be used for this environment
!>        it will not be freed
!> \param ierr will return a number different from 0 if there was an error
!> \note
!>      moved here because of circular dependencies
!> \author fawzi
! *****************************************************************************
  SUBROUTINE run_input(input_file_path,output_file_path,ierr,mpi_comm)
    CHARACTER(len=*), INTENT(in)             :: input_file_path, &
                                                output_file_path
    INTEGER, INTENT(out)                     :: ierr
    INTEGER, INTENT(in), OPTIONAL            :: mpi_comm

    CHARACTER(len=*), PARAMETER :: routineN = 'run_input', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
    LOGICAL                                  :: failure, success
    TYPE(cp_error_type)                      :: error
    TYPE(cp_para_env_type), POINTER          :: para_env

    failure=.FALSE.
    IF (PRESENT(mpi_comm)) THEN
       NULLIFY(para_env)
       CALL cp_para_env_create(para_env, group=mpi_comm, owns_group=.FALSE.,error=error) !XXXXXXXXXXXX uninitiliased error
    ELSE
       para_env => f77_default_para_env
       CALL cp_para_env_retain(para_env,error=error) !XXXXXXXXXXXX uninitiliased error
    END IF
    IF (para_env%mepos==para_env%source) THEN
       IF (output_file_path=="__STD_OUT__") THEN
          unit_nr=6
       ELSE
          CALL open_file(file_name=output_file_path,file_status="UNKNOWN",&
               file_action="WRITE", file_position="APPEND",&
               unit_number=unit_nr)
       END IF
    ELSE
       unit_nr=-1
    END IF
    success=cp2k_run(input_file_path,unit_nr,para_env%group)
    IF (.NOT.success) THEN
       ierr=cp_failure_level
    ELSE
       ierr=0
    END IF
    CALL cp_para_env_release(para_env,error=error) !XXXXXXXXXXXX uninitiliased error
  END SUBROUTINE run_input

END MODULE cp2k_runs
