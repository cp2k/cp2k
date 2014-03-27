!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Tree Monte Carlo entry point, set up, CPU redistribution and 
!>      input reading
!> \par History
!>      11.2012 created [Mandes Schönherr]
!> \author Mandes 
! *****************************************************************************


MODULE tmc_setup
  USE F77_interface,                   ONLY: create_force_env,&
                                             destroy_force_env
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_para_env,                     ONLY: cp_para_env_create
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE environment,                     ONLY: cp2k_get_walltime
  USE global_types,                    ONLY: global_environment_type
  USE header,                          ONLY: tmc_ana_header,&
                                             tmc_header
  USE input_section_types,             ONLY: section_type,&
                                             section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE machine,                         ONLY: default_output_unit,&
                                             m_flush
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_comm_dup,&
                                             mp_comm_split_direct
  USE parallel_rng_types,              ONLY: UNIFORM,&
                                             create_rng_stream,&
                                             delete_rng_stream
  USE physcon,                         ONLY: au2a => angstrom,&
                                             au2bar => bar
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE tmc_analysis,                    ONLY: analysis_init,&
                                             analysis_restart_print,&
                                             analysis_restart_read,&
                                             analyze_file_configurations,&
                                             finalize_tmc_analysis,&
                                             tmc_read_ana_input
  USE tmc_analysis_types,              ONLY: tmc_ana_env_release,&
                                             tmc_ana_list_type
  USE tmc_file_io,                     ONLY: expand_file_name_int
  USE tmc_master,                      ONLY: do_tmc_master
  USE tmc_move_handle,                 ONLY: finalize_mv_types,&
                                             print_move_types,&
                                             read_init_move_types
  USE tmc_stati,                       ONLY: &
       task_type_MC, task_type_ideal_gas, tmc_NMC_worker_out_file_name, &
       tmc_ana_out_file_name, tmc_default_restart_in_file_name, &
       tmc_default_restart_out_file_name, tmc_default_unspecified_name, &
       tmc_energy_worker_out_file_name, tmc_master_out_file_name
  USE tmc_tree_build,                  ONLY: allocate_new_sub_tree_node,&
                                             deallocate_sub_tree_node
  USE tmc_tree_types,                  ONLY: tree_type
  USE tmc_types,                       ONLY: tmc_comp_set_type,&
                                             tmc_env_create,&
                                             tmc_env_release,&
                                             tmc_env_type,&
                                             tmc_master_env_create,&
                                             tmc_master_env_release,&
                                             tmc_worker_env_create,&
                                             tmc_worker_env_release
  USE tmc_worker,                      ONLY: do_tmc_worker,&
                                             get_atom_kinds_and_cell,&
                                             get_initial_conf
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tmc_setup'

  PUBLIC  :: do_tmc, do_analyze_files

CONTAINS

! *****************************************************************************
!> \brief tmc_entry point
!> \param force_env the force environment for the simulation
!> \param globenv the global environment for the simulation
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE do_tmc (input_declaration, root_section, para_env, globenv, error )
    TYPE(section_type), POINTER              :: input_declaration
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'do_tmc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bcast_output_unit, handle, i, &
                                                ierr, output_unit
    LOGICAL                                  :: failure, success
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: init_rng_seed
    TYPE(cp_error_type)                      :: error_sub
    TYPE(cp_logger_type), POINTER            :: logger, logger_sub
    TYPE(section_vals_type), POINTER         :: tmc_ana_section
    TYPE(tmc_ana_list_type), DIMENSION(:), &
      POINTER                                :: tmc_ana_env_list
    TYPE(tmc_env_type), POINTER              :: tmc_env

    failure = .FALSE.
    ! start the timing
    CALL timeset(routineN,handle)

    NULLIFY(logger, logger_sub, tmc_env, tmc_ana_env_list)
    logger=>cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    ! write header, on the 'rank 0' of the global communicator
    IF (output_unit > 0) THEN
       CALL tmc_header(output_unit)
       CALL m_flush(output_unit)
    ENDIF
    ! ugly, we need to know the output unit on source, everywhere, in particular
    ! the tmc master
    IF (output_unit .NE. default_output_unit .AND. output_unit .GT.0) THEN
      WRITE (UNIT=output_unit,FMT="(/,T2,A)") REPEAT("-",79)
      WRITE (UNIT=output_unit,FMT="(/,T2,A)") "The TMC output files are:"
      WRITE (UNIT=output_unit,FMT="(/,T2,A)") "tmc_master.out           the TMC master"
      WRITE (UNIT=output_unit,FMT="(/,T2,A)") "tmc_E_worker_*.out       the worker outputs (energy calculations etc.)"
      WRITE (UNIT=output_unit,FMT="(/,T2,A)") "tmc_ana.out              the analysis output"
      WRITE (UNIT=output_unit,FMT="(/,T2,A)") REPEAT("-",79)
    END IF
    bcast_output_unit=output_unit
    CALL mp_bcast(bcast_output_unit,para_env%source,para_env%group)

    ! create tmc_env
    CALL tmc_env_create(tmc_env, error)
    CALL tmc_preread_input(root_section,tmc_env,error=error)
    CALL tmc_redistributing_cores(tmc_env%tmc_comp_set, para_env, &
                                  ana_on_the_fly=tmc_env%tmc_comp_set%ana_on_the_fly, &
                                  success=success, error=error)

    IF(success) THEN
      ! initialize master and worker environment
      IF(tmc_env%tmc_comp_set%group_nr.EQ.0) THEN
        CALL tmc_master_env_create(tmc_env, error) ! create master env
      ELSE IF(tmc_env%tmc_comp_set%group_nr.NE.0) THEN
        CALL tmc_worker_env_create(tmc_env, error) ! create worker env
      END IF
  
      CALL tmc_read_input(root_section,tmc_env,error=error)
      !CALL init_move_types(tmc_params=tmc_env%params, error=error)
  
      ! init random number generator: use determistic random numbers
      IF(tmc_env%tmc_comp_set%group_nr.EQ.0.AND.&
         tmc_env%m_env%rnd_init.GT.0) THEN
        ALLOCATE(init_rng_seed(3,2))
        init_rng_seed(:,:) = &
            RESHAPE( (/ tmc_env%m_env%rnd_init*42.0_dp, &
                        tmc_env%m_env%rnd_init*54.0_dp, &
                        tmc_env%m_env%rnd_init*63.0_dp, &
                        tmc_env%m_env%rnd_init*98.0_dp, &
                        tmc_env%m_env%rnd_init*10.0_dp, &
                        tmc_env%m_env%rnd_init*2.0_dp/),&
                     (/ 3, 2 /) )
        CALL create_rng_stream(rng_stream=tmc_env%rng_stream, &
                               name="TMC_deterministic_rng_stream",& 
                               seed=init_rng_seed(:,:), &
                               distribution_type=UNIFORM, error=error)
        DEALLOCATE(init_rng_seed)
      ELSE
        CALL create_rng_stream(rng_stream=tmc_env%rng_stream,&
                               name="TMC_rng_stream",&
                               distribution_type=UNIFORM,error=error)
      END IF


      ! start running master and worker routines
      ! the master
      IF(tmc_env%tmc_comp_set%group_nr.EQ.0) THEN
        !TODO get the correct usage of creating and handling the logger...
        CALL cp_logger_create(logger_sub,para_env=tmc_env%tmc_comp_set%para_env_m_only,&
                default_global_unit_nr=default_output_unit, close_global_unit_on_dealloc=.FALSE.)
        CALL cp_logger_set(logger_sub,local_filename="tmc_main")
        CALL cp_error_init(error_sub, stop_level=cp_failure_level, logger=logger_sub)
     
        ! if we're doing output to the screen, keep it there, else this master
        ! opens a file (not that two different ranks are writing to the
        ! default_output_unit, we leave it up to mpirun or so to merge stuff
        IF (bcast_output_unit == default_output_unit) THEN
          tmc_env%m_env%io_unit = default_output_unit
        ELSE
          CALL open_file(file_name=tmc_master_out_file_name,file_status="UNKNOWN",&
               file_action="WRITE", file_position="APPEND",&
               unit_number=tmc_env%m_env%io_unit)
          CALL tmc_header(tmc_env%m_env%io_unit)
        ENDIF
        ! print the intresting parameters and starting values
        CALL tmc_print_params(tmc_env, error_sub)
        CALL print_move_types(init=.TRUE., file_io=tmc_env%m_env%io_unit, &
                              tmc_params=tmc_env%params, error=error_sub)
        CALL do_tmc_master(tmc_env=tmc_env, globenv=globenv, error=error_sub) ! start the master routine

        IF (bcast_output_unit .NE. tmc_env%m_env%io_unit) THEN
          CALL close_file(unit_number=tmc_env%m_env%io_unit)
        ENDIF
  
        CALL cp_error_dealloc_ref(error_sub,error=error)
        CALL cp_logger_release(logger_sub)
  
      ! the worker groups
      ELSE IF(tmc_env%tmc_comp_set%group_nr.GT.0) THEN
        NULLIFY(logger_sub)
        ! create separate logger and error handler for each worker
        CALL cp_logger_create(logger_sub,para_env=tmc_env%tmc_comp_set%para_env_sub_group,&
                default_global_unit_nr=default_output_unit, close_global_unit_on_dealloc=.FALSE.)
        CALL cp_logger_set(logger_sub,local_filename="tmc_localLog")
        CALL cp_error_init(error_sub, stop_level=cp_failure_level, logger=logger_sub)
        tmc_env%w_env%io_unit = default_output_unit
  
        ! energy worker
        IF(tmc_env%tmc_comp_set%group_nr.LE.tmc_env%tmc_comp_set%group_ener_nr) THEN
          CALL create_force_env(new_env_id=tmc_env%w_env%env_id_ener, &
                  input_declaration=input_declaration,&
                  input_path=tmc_env%params%energy_inp_file, &
                  mpi_comm=tmc_env%tmc_comp_set%para_env_sub_group%group, &
                  output_path=TRIM(expand_file_name_int(file_name=tmc_energy_worker_out_file_name, &
                                   ivalue=tmc_env%tmc_comp_set%group_nr, error=error)), & 
                  ierr=ierr)
          CALL cp_assert(ierr.EQ.0, cp_failure_level,cp_assertion_failed, routineP,&
                         "creating force env result in error "//cp_to_string(ierr),&
                         error)
        END IF
        ! worker for configurational change
        IF(tmc_env%params%NMC_inp_file.NE."" .AND. &
           (tmc_env%tmc_comp_set%group_cc_nr.EQ.0 .OR. &
           tmc_env%tmc_comp_set%group_nr.GT.tmc_env%tmc_comp_set%group_ener_nr))THEN
          CALL create_force_env(new_env_id=tmc_env%w_env%env_id_approx, &
                  input_declaration=input_declaration, &
                  input_path=tmc_env%params%NMC_inp_file, &
                  mpi_comm=tmc_env%tmc_comp_set%para_env_sub_group%group, &
                  output_path=TRIM(expand_file_name_int(file_name=tmc_NMC_worker_out_file_name, &
                                   ivalue=tmc_env%tmc_comp_set%group_nr, error=error)), &
                  ierr=ierr)
          CALL cp_assert(ierr.EQ.0, cp_failure_level,cp_assertion_failed, routineP,&
                         "creating approx force env result in error "//cp_to_string(ierr),&
                         error)
        END IF
        CALL do_tmc_worker(tmc_env=tmc_env, error=error_sub) ! start the worker routine
  
        IF(tmc_env%w_env%env_id_ener.GT.0) &
          CALL destroy_force_env(tmc_env%w_env%env_id_ener, ierr)
        IF(tmc_env%w_env%env_id_approx.GT.0) &
          CALL destroy_force_env(tmc_env%w_env%env_id_approx, ierr)
  
        CALL cp_error_dealloc_ref(error_sub,error=error)
        CALL cp_logger_release(logger_sub)

      ! the analysis group
      ELSE IF(ASSOCIATED(tmc_env%tmc_comp_set%para_env_m_ana)) THEN
        ! unused worker groups can do analysis
        NULLIFY(logger_sub)
        ! create separate logger and error handler for each worker
        CALL cp_logger_create(logger_sub,para_env=tmc_env%tmc_comp_set%para_env_m_ana,&
                default_global_unit_nr=default_output_unit, close_global_unit_on_dealloc=.FALSE.)
        tmc_env%w_env%io_unit = default_output_unit
        CALL cp_logger_set(logger_sub,local_filename="tmc_ana_localLog")
        CALL cp_error_init(error_sub, stop_level=cp_failure_level, logger=logger_sub)
        ! if we're doing output to the screen, keep it there, else this master
        ! opens a file (not that two different ranks are writing to the
        ! default_output_unit, we leave it up to mpirun or so to merge stuff
        IF (bcast_output_unit == default_output_unit) THEN
          output_unit = default_output_unit
        ELSE
          CALL open_file(file_name=tmc_ana_out_file_name,file_status="UNKNOWN",&
               file_action="WRITE", file_position="APPEND",&
               unit_number=output_unit)
          CALL tmc_ana_header(output_unit)
        ENDIF

        ALLOCATE(tmc_ana_env_list(tmc_env%params%nr_temp))
        tmc_ana_section => section_vals_get_subs_vals(root_section,"MOTION%TMC%TMC_ANALYSIS",error=error)
        DO i=1, tmc_env%params%nr_temp
          CALL tmc_read_ana_input(tmc_ana_section, tmc_ana_env_list(i)%temp, error)
          tmc_ana_env_list(i)%temp%io_unit = output_unit
        END DO
        CALL do_tmc_worker(tmc_env=tmc_env, ana_list=tmc_ana_env_list, error=error_sub) ! start the worker routine for analysis
        DO i=1, tmc_env%params%nr_temp
          CALL tmc_ana_env_release(tmc_ana_env_list(i)%temp, error)
        END DO
        IF (bcast_output_unit .NE. output_unit) THEN
          CALL close_file(unit_number=tmc_env%m_env%io_unit)
        ENDIF
        CALL cp_error_dealloc_ref(error_sub,error=error)
        CALL cp_logger_release(logger_sub)
      END IF ! unused worker groups have nothing to do
  
      ! delete the random numbers
      CPPrecondition(ASSOCIATED(tmc_env%rng_stream),cp_failure_level,routineP,error,failure)
      CALL delete_rng_stream(tmc_env%rng_stream,error=error)

      ! deallocate the move types
      CALL finalize_mv_types(tmc_env%params, error)
 
      ! finalize master and worker environment
      IF(tmc_env%tmc_comp_set%group_nr.EQ.0) THEN
        CALL tmc_master_env_release(tmc_env, error) ! release master env
      ELSE IF(tmc_env%tmc_comp_set%group_nr.NE.0) THEN
        CALL tmc_worker_env_release(tmc_env, error) ! release worker env
      END IF ! unused worker groups have nothing to do

    ELSE
      IF(tmc_env%params%print_test_output) THEN
        WRITE(output_unit,*)"TMC|NOTenoughProcessorsX= -999"
        WRITE(output_unit,*)"TMC|NOTcalculatedTotal energy: -999"
      END IF
    END IF
    ! finalize / deallocate everthing
    CALL tmc_env_release(tmc_env, error)

    ! end the timing
    CALL timestop(handle)

  END SUBROUTINE do_tmc

! *****************************************************************************
!> \brief analyze TMC trajectory files
!> \param 
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 03.2013
! *****************************************************************************
  SUBROUTINE do_analyze_files(input_declaration, root_section, para_env, globenv, error)
    TYPE(section_type), POINTER              :: input_declaration
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'do_analyze_files', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: comm, dir_ind, handle, &
                                                my_mpi_world, nr_dim, &
                                                output_unit, temp
    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: error_sub
    TYPE(cp_logger_type), POINTER            :: logger, logger_sub
    TYPE(tmc_ana_list_type), DIMENSION(:), &
      POINTER                                :: ana_list
    TYPE(tmc_env_type), POINTER              :: tmc_env
    TYPE(tree_type), POINTER                 :: elem

    NULLIFY(ana_list, tmc_env, elem, logger, logger_sub)
    failure = .FALSE.

    ! start the timing
    CALL timeset(routineN,handle)

    ! create a TMC environment (also to have a params environment)
    CALL tmc_env_create(tmc_env, error)
    ! duplicate communicator 
    CALL mp_comm_dup(para_env%group, my_mpi_world)
    ! -- spiltting communicators
    CALL mp_comm_split_direct(my_mpi_world,comm,para_env%mepos, 0)
    CALL cp_para_env_create(para_env=tmc_env%tmc_comp_set%para_env_m_ana,&
                                group=comm,error=error)
    CALL cp_assert(para_env%num_pe.EQ.1, &
               cp_warning_level,cp_assertion_failed, routineP,&
               "just one out of "//cp_to_string(para_env%num_pe)//"cores is used ",&
               error)
    ! distribute work to availuble cores
    IF(para_env%mepos.EQ.0) THEN
      !TODO get the correct usage of creating and handling the logger...
      logger=>cp_error_get_logger(error)
      output_unit= cp_logger_get_default_io_unit(logger)
      CPPostconditionNoFail(output_unit.GT.0,cp_failure_level,routineP,error)

      CALL cp_logger_create(logger_sub,para_env=tmc_env%tmc_comp_set%para_env_m_ana,&
             default_global_unit_nr=default_output_unit, close_global_unit_on_dealloc=.FALSE.)
      CALL cp_logger_set(logger_sub,local_filename="tmc_ana")
      CALL cp_error_init(error_sub, stop_level=cp_failure_level, logger=logger_sub)
      ! write the header
      CALL tmc_ana_header(output_unit)

      ! read the input and create the ana environments for each temp
      CALL tmc_read_ana_files_input(input_declaration=input_declaration,&
                                    input=root_section,ana_list=ana_list,&
                                    elem=elem, tmc_env=tmc_env,error=error)
      nr_dim = SIZE(elem%pos)
      ! we need a new tree element with all neccessay arrays, (e.g. dipoles could not be allocated already)
      CALL deallocate_sub_tree_node(tmc_params=tmc_env%params, &
               tree_elem=elem, error=error)
      CPPostconditionNoFail(SIZE(ana_list).GT.0,cp_failure_level,routineP,error)

      ! print initial test output (for single core tests, where no data is produced)
      IF(tmc_env%params%print_test_output) THEN
        WRITE(output_unit,*)"TMC|ANAtestOutputInitX= -999"
      END IF

      ! do the analysis
      DO temp=1, SIZE(ana_list)
        ! initialize the structures
        ana_list(temp)%temp%io_unit = output_unit
        CALL analysis_init(ana_env=ana_list(temp)%temp, nr_dim=nr_dim, error=error)
        ! to allocate the dipole array in tree elements
        IF(ana_list(temp)%temp%costum_dip_file_name.NE.&
           tmc_default_unspecified_name)&
          tmc_env%params%print_dipole = .TRUE.

        IF(.NOT.ASSOCIATED(elem)) &
          CALL allocate_new_sub_tree_node(tmc_params=tmc_env%params, &
                 next_el=elem, nr_dim=nr_dim, &
                 error=error)
        CALL analysis_restart_read(ana_env=ana_list(temp)%temp, &
                                   elem=elem, error=error)
        CALL cp_assert(ASSOCIATED(elem).OR.&
               ASSOCIATED(ana_list(temp)%temp%last_elem), &
               cp_failure_level,cp_assertion_failed, routineP,&
               "uncorrect initialization of the initial configuration ",&
               error)
        ! do for all directories
        DO dir_ind=1, SIZE(ana_list(temp)%temp%dirs)
          WRITE(output_unit,FMT='(T2,A,"| ",A,T41,A40)')"TMC_ANA",&
            "read directory",TRIM(ana_list(temp)%temp%dirs(dir_ind))
          CALL analyze_file_configurations(&
                 start_id=ana_list(temp)%temp%from_elem, &
                 end_id=ana_list(temp)%temp%to_elem, &
                 dir_ind=dir_ind, &
                 ana_env=ana_list(temp)%temp,&
                 tmc_params=tmc_env%params, error=error)
          ! remove the last saved element to start with a new file
          !  there is no weight for this element
          IF(dir_ind.LT.SIZE(ana_list(temp)%temp%dirs) .AND. &
             ASSOCIATED(ana_list(temp)%temp%last_elem)) &
            CALL deallocate_sub_tree_node(tmc_params=tmc_env%params, &
                   tree_elem=ana_list(temp)%temp%last_elem, error=error)
          IF(ASSOCIATED(ana_list(temp)%temp%last_elem)) &
            ana_list(temp)%temp%conf_offset = ana_list(temp)%temp%conf_offset&
              + ana_list(temp)%temp%last_elem%nr
        END DO
        CALL finalize_tmc_analysis(ana_env=ana_list(temp)%temp, error=error)
        ! write analysis restart file
        !  if there is somthing to write
        ! shifts the last element to actual element 
        IF(ASSOCIATED(ana_list(temp)%temp%last_elem)) & 
          CALL analysis_restart_print(ana_env=ana_list(temp)%temp, &
                                      error=error)
        IF(ASSOCIATED(ana_list(temp)%temp%last_elem)) &
          CALL deallocate_sub_tree_node(tmc_params=tmc_env%params, &
                 tree_elem=ana_list(temp)%temp%last_elem, error=error)
        IF(ASSOCIATED(elem)) &
          CALL deallocate_sub_tree_node(tmc_params=tmc_env%params, &
                 tree_elem=elem, error=error)

        CALL tmc_ana_env_release(ana_list(temp)%temp, error)
      END DO

      DEALLOCATE(ana_list)
  
      ! finalize / deallocate everthing
      CALL cp_error_dealloc_ref(error_sub,error=error)
      CALL cp_logger_release(logger_sub)
    END IF
    CALL tmc_env_release(tmc_env, error)

    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE do_analyze_files

! *****************************************************************************
!> \brief creates a new para environment for tmc analysis for each temperature
!> \param input global environment
!> \param tmc_ana TMC analysis environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 03.2013
! *****************************************************************************
  SUBROUTINE tmc_read_ana_files_input(input_declaration, input, ana_list, elem, tmc_env, error)
    TYPE(section_type), POINTER              :: input_declaration
    TYPE(section_vals_type), POINTER         :: input
    TYPE(tmc_ana_list_type), DIMENSION(:), &
      POINTER                                :: ana_list
    TYPE(tree_type), POINTER                 :: elem
    TYPE(tmc_env_type), POINTER              :: tmc_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tmc_read_ana_files_input', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length), &
      DIMENSION(:), POINTER                  :: directories
    INTEGER                                  :: env_id, ierr, nr_temp, t_act
    LOGICAL                                  :: failure, flag
    REAL(KIND=dp)                            :: tmax, tmin
    REAL(KIND=dp), DIMENSION(:), POINTER     :: inp_Temp, Temps
    TYPE(section_vals_type), POINTER         :: tmc_section

    NULLIFY(tmc_section, inp_Temp, Temps)
    failure = .FALSE.
    CPPrecondition(ASSOCIATED(input),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(ana_list),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_env),cp_failure_level,routineP,error,failure)

    ! first global TMC stuff
    tmc_section => section_vals_get_subs_vals(input,"MOTION%TMC",error=error)
    CALL section_vals_val_get(tmc_section,"PRINT_TEST_OUTPUT",l_val=tmc_env%params%print_test_output,error=error)
    ! TMC analysis stuff
    tmc_section => section_vals_get_subs_vals(input,"MOTION%TMC%TMC_ANALYSIS_FILES",error=error)
    CALL section_vals_get(tmc_section,explicit=flag,error=error)
    CPPostconditionNoFail(flag,cp_failure_level,routineP,error)

    CALL section_vals_val_get(tmc_section,"FORCE_ENV_FILE",&
             c_val=tmc_env%params%energy_inp_file,error=error)

    CALL section_vals_val_get(tmc_section,"NR_TEMPERATURE",i_val=nr_temp,error=error)

    CALL section_vals_val_get(tmc_section,"TEMPERATURE",r_vals=inp_Temp,error=error)
    CALL cp_assert(.NOT.((nr_temp .GT. 1).AND.(SIZE(inp_Temp).NE.2)),&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "specify each temperature, skip keyword NR_TEMPERATURE",&
                     error, failure=flag, only_ionode=.TRUE.)
    IF(nr_temp .EQ. 1) THEN
      nr_temp = SIZE(inp_Temp)
      ALLOCATE(Temps(nr_temp))
      Temps(:) = inp_Temp(:)
    ELSE
      tmin   = inp_Temp(1)
      tmax   = inp_Temp(2)
      ALLOCATE(Temps(nr_temp))
      Temps(1) = tmin
      DO t_act = 2, SIZE(Temps)
        Temps(t_act) = Temps(t_act-1) + (tmax-tmin)/(SIZE(Temps)-1.0_dp)
      END DO
      CALL cp_assert(ALL(Temps.GE.0.0_dp),&
                   cp_failure_level,cp_assertion_failed,routineP,&
                   "The temperatures are negative. Should be specified using "//&
                   "TEMPERATURE {T_min} {T_max} and NR_TEMPERATURE {#temperatures}",&
                   error, failure=flag, only_ionode=.TRUE.)
    END IF

    ! get multiple directories
    CALL section_vals_val_get(tmc_section,"DIRECTORIES",c_vals=directories,error=error)

    ! get init configuration (for sizes)
    CALL create_force_env(new_env_id=env_id, &
           input_declaration=input_declaration,&
           input_path=tmc_env%params%energy_inp_file, &
           mpi_comm=tmc_env%tmc_comp_set%para_env_m_ana%group, &
           output_path="tmc_ana.out", ierr=ierr)
    CALL get_initial_conf(tmc_params=tmc_env%params, init_conf=elem, &
                          env_id=env_id, error=error)
    CALL get_atom_kinds_and_cell(env_id=env_id, atoms=tmc_env%params%atoms, &
                                 cell=tmc_env%params%cell, error=error)
    CALL destroy_force_env(env_id, ierr)
    
    ALLOCATE(ana_list(SIZE(Temps)))
    DO t_act=1, SIZE(Temps)
      ana_list(t_act)%temp => NULL()
      CALL tmc_read_ana_input(tmc_section, ana_list(t_act)%temp, error)
      ana_list(t_act)%temp%temperature = Temps(t_act)
      ALLOCATE(ana_list(t_act)%temp%dirs(SIZE(directories)))
      ana_list(t_act)%temp%dirs(:)     =  directories(:)
      ana_list(t_act)%temp%cell        => tmc_env%params%cell
      ana_list(t_act)%temp%atoms       => tmc_env%params%atoms
      ana_list(t_act)%temp%print_test_output = tmc_env%params%print_test_output

      CALL section_vals_val_get(tmc_section,"POSITION_FILE",&
             c_val=ana_list(t_act)%temp%costum_pos_file_name,error=error)
      CALL section_vals_val_get(tmc_section,"DIPOLE_FILE",&
             c_val=ana_list(t_act)%temp%costum_dip_file_name,error=error)
      CALL section_vals_val_get(tmc_section,"CELL_FILE",&
             c_val=ana_list(t_act)%temp%costum_cell_file_name,error=error)
      CALL section_vals_val_get(tmc_section,"START_ELEM",i_val=ana_list(t_act)%temp%from_elem,error=error)
      CALL section_vals_val_get(tmc_section,"END_ELEM",i_val=ana_list(t_act)%temp%to_elem,error=error)
    END DO
    DEALLOCATE(Temps)
  END SUBROUTINE tmc_read_ana_files_input

! *****************************************************************************
!> \brief read the variables for distributing cores
!> \param tmc_section points to the tmc subsection in the input file
!> \param tmc_env structure for storing all the tmc parameters
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE tmc_preread_input(input,tmc_env, error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(tmc_env_type), POINTER              :: tmc_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tmc_preread_input', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: c_tmp
    INTEGER                                  :: itmp
    LOGICAL                                  :: explicit_key, failure, flag
    REAL(KIND=dp)                            :: tmax, tmin
    REAL(KIND=dp), DIMENSION(:), POINTER     :: inp_Temp
    TYPE(section_vals_type), POINTER         :: tmc_section

    failure = .FALSE.
    NULLIFY(tmc_section, inp_Temp)

    CPPrecondition(ASSOCIATED(input),cp_failure_level,routineP,error,failure)

    tmc_env%tmc_comp_set%ana_on_the_fly = 0
    tmc_section => section_vals_get_subs_vals(input,"MOTION%TMC%TMC_ANALYSIS",error=error)
    CALL section_vals_get(tmc_section,explicit=flag,error=error)
    IF(flag) THEN
      tmc_env%tmc_comp_set%ana_on_the_fly = 1
    END IF

    tmc_section => section_vals_get_subs_vals(input,"MOTION%TMC",error=error)
    CALL section_vals_get(tmc_section,explicit=flag,error=error)
    CPPostconditionNoFail(flag,cp_failure_level,routineP,error)

    CALL section_vals_val_get(tmc_section,"PRINT_TEST_OUTPUT",l_val=tmc_env%params%print_test_output,error=error)

    CPPrecondition(ASSOCIATED(tmc_env%tmc_comp_set),cp_failure_level,routineP,error,failure)
    ! read the parameters for the computational setup
    CALL section_vals_val_get(tmc_section,"GROUP_ENERGY_SIZE",i_val=tmc_env%tmc_comp_set%group_ener_size,error=error)
    CALL section_vals_val_get(tmc_section,"GROUP_ENERGY_NR",i_val=tmc_env%tmc_comp_set%group_ener_nr,error=error)
    CALL section_vals_val_get(tmc_section,"GROUP_CC_SIZE",i_val=tmc_env%tmc_comp_set%group_cc_size,error=error)
    CALL section_vals_val_get(tmc_section,"GROUP_ANLYSIS_NR",i_val=itmp, error=error)
    IF(tmc_env%tmc_comp_set%ana_on_the_fly.GT.0) &
      tmc_env%tmc_comp_set%ana_on_the_fly = itmp
    CALL cp_assert(tmc_env%tmc_comp_set%ana_on_the_fly.LE.1,&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "analysing on the fly is up to now not supported for multiple cores. "//&
                     " Restart file witing for this case and temperature "//&
                     "distribution has to be solved.!.",&
                     error, failure=flag, only_ionode=.TRUE.)
    CALL section_vals_val_get(tmc_section,"RESULT_LIST_IN_MEMORY",l_val=tmc_env%params%USE_REDUCED_TREE,error=error)
    ! swap the variable, because of oposit meaning
    tmc_env%params%USE_REDUCED_TREE = .NOT.tmc_env%params%USE_REDUCED_TREE
    CALL section_vals_val_get(tmc_section,"NR_TEMPERATURE",i_val=tmc_env%params%nr_temp,error=error)

    ! stuff everyone needs to know
    CALL section_vals_val_get(tmc_section,"NMC_MOVES%NMC_FILE_NAME",c_val=tmc_env%params%NMC_inp_file,error=error)
    IF(tmc_env%params%NMC_inp_file.EQ.tmc_default_unspecified_name)THEN
    ! file name keyword without file name
      CALL cp_assert(.FALSE.,&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "no or a valid NMC input file has to be specified ",&
                     error, failure=flag, only_ionode=.TRUE.)
    ELSE IF(tmc_env%params%NMC_inp_file.EQ."") THEN
    ! no keyword
      CALL cp_assert(tmc_env%tmc_comp_set%group_cc_size.LE.0,&
                     cp_warning_level,cp_assertion_failed,routineP,&
                     "The configurational groups are deactivated, "//&
                     "because no approximated energy input is specified.",&
                     error, failure=flag, only_ionode=.TRUE.)
      tmc_env%tmc_comp_set%group_cc_size = 0
    ELSE
    ! check file existance
      INQUIRE (FILE=TRIM(tmc_env%params%NMC_inp_file),EXIST=flag,IOSTAT=itmp)
      CALL cp_assert(flag .AND. itmp.EQ.0,&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "a valid NMC input file has to be specified ",&
                     error, failure=flag, only_ionode=.TRUE.)
    END IF

    CALL section_vals_val_get(tmc_section,"TEMPERATURE",r_vals=inp_Temp,error=error)
    CALL cp_assert(.NOT.((tmc_env%params%nr_temp .GT. 1).AND.(SIZE(inp_Temp).NE.2)),&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "specify each temperature, skip keyword NR_TEMPERATURE",&
                     error, failure=flag, only_ionode=.TRUE.)
    IF(tmc_env%params%nr_temp .EQ. 1) THEN
      tmc_env%params%nr_temp = SIZE(inp_Temp)
      ALLOCATE(tmc_env%params%Temp(tmc_env%params%nr_temp))
      tmc_env%params%Temp(:) = inp_Temp(:)
    ELSE    
      tmin   = inp_Temp(1)
      tmax   = inp_Temp(2)
      ALLOCATE(tmc_env%params%Temp(tmc_env%params%nr_temp))
      tmc_env%params%Temp(1) = tmin
      DO itmp = 2, SIZE(tmc_env%params%Temp)
        tmc_env%params%Temp(itmp) = tmc_env%params%Temp(itmp-1) + (tmax-tmin)/(SIZE(tmc_env%params%Temp)-1.0_dp)
      END DO
      CALL cp_assert(ALL(tmc_env%params%Temp.GE.0.0_dp),&
                   cp_failure_level,cp_assertion_failed,routineP,&
                   "The temperatures are negative. Should be specified using "//&
                   "TEMPERATURE {T_min} {T_max} and NR_TEMPERATURE {#temperatures}",&
                   error, failure=flag, only_ionode=.TRUE.)
    END IF

    CALL section_vals_val_get(tmc_section,"TASK_TYPE" ,explicit=explicit_key,error=error)
    IF(explicit_key)THEN
      CALL section_vals_val_get(tmc_section,"TASK_TYPE",c_val=c_tmp,error=error)
      SELECT CASE (TRIM(c_tmp))
      CASE (TRIM(tmc_default_unspecified_name))
        tmc_env%params%task_type = task_type_MC
      CASE ("IDEAL_GAS")
        tmc_env%params%task_type = task_type_ideal_gas
      CASE DEFAULT
        CALL cp_assert(.FALSE.,&
               cp_warning_level,cp_assertion_failed,routineP,&
               'unknown TMC task type "'//TRIM(c_tmp)//'" specified. '//&
               " Set to default.",&
               error, failure=flag, only_ionode=.TRUE.)
        tmc_env%params%task_type = task_type_MC
      END SELECT
    END IF

  END SUBROUTINE tmc_preread_input


! *****************************************************************************
!> \brief read the tmc subsection from the input file
!> \param tmc_section points to the tmc subsection in the input file
!> \param tmc_env structure for storing all the tmc parameters
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE tmc_read_input(input,tmc_env, error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(tmc_env_type), POINTER              :: tmc_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tmc_read_input', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: itmp
    LOGICAL                                  :: explicit, failure, flag
    REAL(KIND=dp)                            :: r_tmp
    REAL(KIND=dp), DIMENSION(:), POINTER     :: r_arr_tmp
    TYPE(section_vals_type), POINTER         :: tmc_section

    NULLIFY(tmc_section)

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(input),cp_failure_level,routineP,error,failure)

    tmc_section => section_vals_get_subs_vals(input,"MOTION%TMC",error=error)
    CALL section_vals_get(tmc_section,explicit=flag,error=error)
    CPPostconditionNoFail(flag,cp_failure_level,routineP,error)

    ! only for the master
    IF(tmc_env%tmc_comp_set%group_nr == 0) THEN 
      CPPrecondition(ASSOCIATED(tmc_env%m_env),cp_failure_level,routineP,error,failure)
      ! the walltime input can be done as HH:MM:SS or just in seconds.
      CALL cp2k_get_walltime(section=input, keyword_name="GLOBAL%WALLTIME", &
                             walltime=tmc_env%m_env%walltime, error=error)

      CALL section_vals_val_get(tmc_section,"NUM_MC_ELEM",i_val=tmc_env%m_env%num_MC_elem,error=error)
      CALL section_vals_val_get(tmc_section,"RND_DETERMINISTIC",i_val=tmc_env%m_env%rnd_init,error=error)
      ! restarting
      CALL section_vals_val_get(tmc_section,"RESTART_IN",c_val=tmc_env%m_env%restart_in_file_name,error=error)
      IF(tmc_env%m_env%restart_in_file_name.EQ.tmc_default_unspecified_name) THEN
        tmc_env%m_env%restart_in_file_name = tmc_default_restart_in_file_name
        INQUIRE(FILE=tmc_env%m_env%restart_in_file_name, EXIST=flag)
        IF(.NOT.flag) tmc_env%m_env%restart_in_file_name = ""
      END IF
      CALL section_vals_val_get(tmc_section,"RESTART_OUT",i_val=tmc_env%m_env%restart_out_step,error=error)
      ! restart just at the end (lone keyword)
      IF(tmc_env%m_env%restart_out_step.EQ.-9) THEN
        tmc_env%m_env%restart_out_file_name = tmc_default_restart_out_file_name
        tmc_env%m_env%restart_out_step = HUGE(tmc_env%m_env%restart_out_step)
      END IF
      CALL cp_assert(tmc_env%m_env%restart_out_step.GE.0,&
                   cp_failure_level,cp_assertion_failed,routineP,&
                   "Please specify a valid value for the frequency "//&
                   "to write restart files (RESTART_OUT #). "//&
                   "# > 0 to define the amount of Markov chain elements in between, "//&
                   "or 0 to deactivate the restart file writing. "//&
                   "Lonely keyword writes restart file only at the end of the run.",&
                   error, failure=flag, only_ionode=.TRUE.)
 
      CALL section_vals_val_get(tmc_section,"INFO_OUT_STEP_SIZE",i_val=tmc_env%m_env%info_out_step_size,error=error)
      CALL section_vals_val_get(tmc_section,"DOT_TREE",c_val=tmc_env%params%dot_file_name,error=error)
      CALL section_vals_val_get(tmc_section,"ALL_CONF_FILE_NAME",c_val=tmc_env%params%all_conf_file_name,error=error)
      IF(tmc_env%params%dot_file_name.NE."") tmc_env%params%DRAW_TREE = .TRUE.

    ! everything for the worker group
    ELSE IF(tmc_env%tmc_comp_set%group_nr .NE. 0) THEN 
      CPPrecondition(ASSOCIATED(tmc_env%w_env),cp_failure_level,routineP,error,failure)
    END IF

    ! stuff everyone needs to know

    ! the NMC_FILE_NAME is already read in tmc_preread_input
    CALL section_vals_val_get(tmc_section,"ENERGY_FILE_NAME",c_val=tmc_env%params%energy_inp_file,error=error)
    ! file name keyword without file name
    CALL cp_assert(tmc_env%params%energy_inp_file.NE."",&
                   cp_failure_level,cp_assertion_failed,routineP,&
                   "a valid exact energy input file has to be specified ",&
                   error, failure=flag, only_ionode=.TRUE.)
    ! check file existance
    INQUIRE (FILE=TRIM(tmc_env%params%energy_inp_file),EXIST=flag,IOSTAT=itmp)
    CALL cp_assert(flag.AND.itmp.EQ.0,&
                   cp_failure_level,cp_assertion_failed,routineP,&
                   "a valid exact energy input file has to be specified, "//&
                   TRIM(tmc_env%params%energy_inp_file)//" does not exist.",&
                   error, failure=flag, only_ionode=.TRUE.)

    CALL section_vals_val_get(tmc_section,"NUM_MV_ELEM_IN_CELL",i_val=tmc_env%params%nr_elem_mv,error=error)    

    CALL section_vals_val_get(tmc_section,"VOLUME_ISOTROPIC",l_val=tmc_env%params%v_isotropic,error=error)
    CALL section_vals_val_get(tmc_section,"PRESSURE",r_val=tmc_env%params%pressure,error=error)
    tmc_env%params%pressure = tmc_env%params%pressure/au2bar
    CALL section_vals_val_get(tmc_section,"MOVE_CENTER_OF_MASS",l_val=tmc_env%params%mv_cen_of_mass,error=error)

    CALL section_vals_val_get(tmc_section,"SUB_BOX",r_vals=r_arr_tmp,error=error)
    IF(SIZE(r_arr_tmp).GT.1) THEN
      CALL cp_assert(SIZE(r_arr_tmp).EQ.tmc_env%params%dim_per_elem,&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "The entered sub box sizes does not fit in number of dimensions.",&
                     error, failure=flag, only_ionode=.TRUE.)
      CALL cp_assert(ALL(r_arr_tmp.GT.0.0_dp),&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "The entered sub box lenghts should be greater than 0.",&
                     error, failure=flag, only_ionode=.TRUE.)
      DO itmp=1, SIZE(tmc_env%params%sub_box_size)
        tmc_env%params%sub_box_size(itmp) = r_arr_tmp(itmp)/au2a
      END DO
    ELSE IF(r_arr_tmp(1).GT.0.0_dp) THEN
      r_tmp = r_arr_tmp(1)/au2a
      tmc_env%params%sub_box_size(:) = r_tmp
    END IF

    ! read all the distinct moves
    CALL read_init_move_types(tmc_params=tmc_env%params, &
                              tmc_section=tmc_section, error=error)

    CALL section_vals_val_get(tmc_section,"ESIMATE_ACC_PROB",l_val=tmc_env%params%esimate_acc_prob,error=error)
    CALL section_vals_val_get(tmc_section,"SPECULATIVE_CANCELING",l_val=tmc_env%params%SPECULATIVE_CANCELING,error=error)
    CALL section_vals_val_get(tmc_section,"USE_SCF_ENERGY_INFO",l_val=tmc_env%params%use_scf_energy_info,error=error)
    ! printing
    CALL section_vals_val_get(tmc_section,"PRINT_ONLY_ACC",l_val=tmc_env%params%print_only_diff_conf,error=error)
    CALL section_vals_val_get(tmc_section,"PRINT_COORDS",l_val=tmc_env%params%print_trajectory,error=error)
    CALL section_vals_val_get(tmc_section,"PRINT_DIPOLE" ,explicit=explicit,error=error)
    IF(explicit) &
      CALL section_vals_val_get(tmc_section,"PRINT_DIPOLE",l_val=tmc_env%params%print_dipole,error=error)
    CALL section_vals_val_get(tmc_section,"PRINT_FORCES" ,explicit=explicit,error=error)
    IF(explicit) &
      CALL section_vals_val_get(tmc_section,"PRINT_FORCES",l_val=tmc_env%params%print_forces,error=error)
    CALL section_vals_val_get(tmc_section,"PRINT_CELL" ,explicit=explicit,error=error)
    IF(explicit) &
      CALL section_vals_val_get(tmc_section,"PRINT_CELL",l_val=tmc_env%params%print_cell,error=error)
    CALL section_vals_val_get(tmc_section,"PRINT_ENERGIES",l_val=tmc_env%params%print_energies,error=error)

  END SUBROUTINE tmc_read_input

! *****************************************************************************
!> \brief creates a new para environment for tmc
!> \param tmc_comp_set structure with parameters for computational setup
!> \param para_env the old parallel environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE tmc_redistributing_cores(tmc_comp_set,para_env, ana_on_the_fly, &
                                      success, error)
    TYPE(tmc_comp_set_type), POINTER         :: tmc_comp_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: ana_on_the_fly
    LOGICAL                                  :: success
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tmc_redistributing_cores', &
      routineP = moduleN//':'//routineN

    INTEGER :: cc_group, cc_group_rank, comm_tmp, master_ana_group, &
      master_ana_rank, master_first_e_worker_g, master_first_e_worker_r, &
      master_worker_group, master_worker_rank, my_mpi_undefined, &
      my_mpi_world, total_used
    LOGICAL                                  :: failure, flag, master

    CPPrecondition(ASSOCIATED(tmc_comp_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)

    ! colors and positions for new communicators
    ! variables for printing
    tmc_comp_set%group_nr = -1
    my_mpi_undefined =  para_env%num_pe + 10000 !HUGE(my_mpi_undefined)! mp_undefined
    master_worker_group   = my_mpi_undefined
    master_worker_rank    = -1
    cc_group              = my_mpi_undefined
    cc_group_rank         = -1
    master_first_e_worker_g = my_mpi_undefined
    master_first_e_worker_r = -1
    master_ana_group = my_mpi_undefined
    master_ana_rank  = -1

    master = .FALSE.
    failure = .FALSE.
    flag = .FALSE.
    success = .TRUE.

    IF(para_env%num_pe.LE.1) THEN
      CALL cp_assert(.FALSE.,&
                     cp_warning_level,cp_assertion_failed,routineP,&
                     "TMC need at least 2 cores (one for master, one for worker)",&
                     error, failure=flag, only_ionode=.TRUE.) 
      success = .FALSE.
    ELSE
      ! check if there are enougth cores available
      CALL cp_assert(tmc_comp_set%group_ener_size*tmc_comp_set%group_ener_nr.LE.(para_env%num_pe-1),&
                     cp_warning_level,cp_assertion_failed,routineP,&
                     "The selected energy group size is too huge. ",&
                     error, failure=flag, only_ionode=.TRUE.)
      IF(flag)THEN
        tmc_comp_set%group_ener_nr = INT((para_env%num_pe-1)/&
          REAL(tmc_comp_set%group_ener_size,KIND=dp))
        CALL cp_assert(tmc_comp_set%group_ener_nr.GE.1,&
                     cp_warning_level,cp_assertion_failed,routineP,&
                     "The selected energy group size is too huge. ",&
                     error, failure=flag, only_ionode=.TRUE.)
        IF(flag) success = .FALSE.
      END IF
  
      ! set the amount of configurational change worker groups
      tmc_comp_set%group_cc_nr = 0 
      IF(tmc_comp_set%group_cc_size.GT.0)THEN
        tmc_comp_set%group_cc_nr = INT((para_env%num_pe-1 -tmc_comp_set%ana_on_the_fly&
            - tmc_comp_set%group_ener_size*tmc_comp_set%group_ener_nr)/&
                                       REAL(tmc_comp_set%group_cc_size,KIND=dp))
  
        CALL cp_assert(tmc_comp_set%group_cc_nr.GE.1,&
                       cp_warning_level,cp_assertion_failed,routineP,&
                       "There are not enougth cores left for creating groups for configurational change.",&
                       error, failure=flag, only_ionode=.TRUE.)
        IF(flag) success = .FALSE.
      END IF
  
      total_used = tmc_comp_set%group_ener_size*tmc_comp_set%group_ener_nr+&
                     tmc_comp_set%group_cc_size*tmc_comp_set%group_cc_nr + &
                     tmc_comp_set%ana_on_the_fly
      CALL cp_assert(para_env%num_pe-1.LE.total_used,&
                     cp_warning_level,cp_assertion_failed,routineP,&
                     cp_to_string(para_env%num_pe-1-total_used)//&
                     " mpi ranks are unused, but can be used for analysis.",&
                     error, failure=flag, only_ionode=.TRUE.)
  
      ! duplicate communicator 
      CALL mp_comm_dup(para_env%group, my_mpi_world)
  
      ! determine the master node
      IF(para_env%mepos==para_env%num_pe-1)THEN
        master = .TRUE.
        master_worker_group     = para_env%num_pe +3 ! belong to master_worker_comm
        master_worker_rank      = 0                  ! rank in m_w_comm
        master_first_e_worker_g = para_env%num_pe +3 ! belong to master_first_energy_worker_comm
        master_first_e_worker_r = 0 
        tmc_comp_set%group_nr   = 0 !para_env%num_pe +3
        master_ana_group        = para_env%num_pe +4
        master_ana_rank         = 0
      ELSE
        ! energy calculation groups
        IF(para_env%mepos .LT. tmc_comp_set%group_ener_size*tmc_comp_set%group_ener_nr)THEN
          tmc_comp_set%group_nr = INT(para_env%mepos/tmc_comp_set%group_ener_size)+1 ! assign to groups
          ! master of worker group
          IF(MODULO(para_env%mepos,tmc_comp_set%group_ener_size).EQ.0) THEN          ! tmc_comp_set%group_nr masters
            master_worker_group = para_env%num_pe +3                                 ! belong to master_worker_comm
            master_worker_rank = tmc_comp_set%group_nr                               ! rank in m_w_comm
            IF(master_worker_rank.EQ.1) THEN
              master_first_e_worker_g = para_env%num_pe +3 ! belong to master_first_energy_worker_comm
              master_first_e_worker_r = 1
            END IF
          END IF
          cc_group = tmc_comp_set%group_nr
          cc_group_rank = para_env%mepos-&
            (tmc_comp_set%group_nr-1)*tmc_comp_set%group_ener_size                   ! rank in worker group
  
        ! configurational change groups
        ELSE IF(para_env%mepos .LT. (tmc_comp_set%group_ener_size*tmc_comp_set%group_ener_nr + &
                                     tmc_comp_set%group_cc_size*tmc_comp_set%group_cc_nr)) THEN
          cc_group_rank = para_env%mepos - tmc_comp_set%group_ener_size*tmc_comp_set%group_ener_nr ! temporary
          tmc_comp_set%group_nr = tmc_comp_set%group_ener_nr+1 + INT(cc_group_rank/tmc_comp_set%group_cc_size) 
          cc_group = tmc_comp_set%group_nr
          ! master of worker group
          IF(MODULO(cc_group_rank,tmc_comp_set%group_cc_size).EQ.0)THEN ! tmc_comp_set%group_nr masters
            master_worker_group = para_env%num_pe +3                    ! belong to master_worker_comm
            master_worker_rank  = tmc_comp_set%group_nr                 ! rank in m_w_comm
          END IF
          !cc_group_rank = cc_group_rank-(tmc_comp_set%group_nr-1)*tmc_comp_set%group_cc_size       ! rank in worker group
          cc_group_rank = MODULO(cc_group_rank, tmc_comp_set%group_cc_size)       ! rank in worker group
        ELSE
          ! not used cores
          ! up to now we use just one core for doing the analysis
          IF(para_env%mepos .EQ. para_env%num_pe-2) THEN
            tmc_comp_set%group_nr = para_env%mepos - (para_env%num_pe-1) ! negative
            CPPostcondition(tmc_comp_set%group_nr.LT.0,cp_failure_level,routineP,error,failure)
            IF(para_env%mepos.GE.para_env%num_pe-1-ana_on_the_fly) THEN
              master_ana_group = para_env%num_pe +4
              master_ana_rank  = -tmc_comp_set%group_nr
            END IF
          END IF
        END IF
      END IF

      IF(success) THEN
        ! -- spiltting communicators
        ! worker intern communication
        CALL mp_comm_split_direct(my_mpi_world,comm_tmp,cc_group, cc_group_rank)
        NULLIFY(tmc_comp_set%para_env_sub_group)
        ! not the unused cores
        IF(cc_group_rank.GE.0) &
          CALL cp_para_env_create(para_env=tmc_comp_set%para_env_sub_group,&
                                  group=comm_tmp,error=error)
  
        ! worker master communication
        CALL mp_comm_split_direct(my_mpi_world,comm_tmp, master_worker_group, master_worker_rank)
        NULLIFY(tmc_comp_set%para_env_m_w)
        ! not the unused cores
        IF(master_worker_rank.GE.0) &
          CALL cp_para_env_create(para_env=tmc_comp_set%para_env_m_w,&
                                  group=comm_tmp,error=error)
  
        ! communicator only for first energy worker master and global master 
        CALL mp_comm_split_direct(my_mpi_world,comm_tmp, master_first_e_worker_g, master_first_e_worker_r)
        NULLIFY(tmc_comp_set%para_env_m_first_w)
        ! not the unused cores
        IF(master_first_e_worker_r.GE.0) &
          CALL cp_para_env_create(para_env=tmc_comp_set%para_env_m_first_w,&
                                  group=comm_tmp,error=error)
   
        ! communicator only for analysis worker and global master
        CALL mp_comm_split_direct(my_mpi_world,comm_tmp, master_ana_group, master_ana_rank)
        NULLIFY(tmc_comp_set%para_env_m_ana)
        IF(master_ana_rank.GE.0) &
          CALL cp_para_env_create(para_env=tmc_comp_set%para_env_m_ana,&
                                  group=comm_tmp,error=error)

        ! communicator for master only to handle external control
        master_ana_group = my_mpi_undefined
        master_ana_rank  = -1
        IF(master) THEN
          master_ana_group = 1
          master_ana_rank  = 1
        END IF
        CALL mp_comm_split_direct(my_mpi_world,comm_tmp, master_ana_group, master_ana_rank)
        NULLIFY(tmc_comp_set%para_env_m_only)
        IF(master_ana_rank.GE.0) &
          CALL cp_para_env_create(para_env=tmc_comp_set%para_env_m_only,&
                                  group=comm_tmp,error=error)
      END IF
    END IF
  END SUBROUTINE tmc_redistributing_cores

! *****************************************************************************
!> \brief prints the most important parameters used for TMC
!> \param tmc_enc structure with parameters for TMC
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE tmc_print_params(tmc_env, error)
    TYPE(tmc_env_type), POINTER              :: tmc_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: fmt_my = '(T2,A,"| ",A,T41,A40)', &
      plabel = "TMC", routineN = 'tmc_print_params', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=80)                        :: c_tmp, fmt_tmp
    INTEGER                                  :: file_nr, i
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(tmc_env),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(tmc_env%tmc_comp_set),cp_failure_level,routineP,error,failure)
    ! only the master prints out
    IF(tmc_env%tmc_comp_set%group_nr == 0) THEN
      file_nr = tmc_env%m_env%io_unit 
      CPPostcondition(ASSOCIATED(tmc_env%tmc_comp_set%para_env_m_w),cp_failure_level,routineP,error,failure)
      CPPostcondition(ASSOCIATED(tmc_env%m_env),cp_failure_level,routineP,error,failure)

      CALL m_flush(file_nr)
      WRITE ( file_nr, * )

      WRITE (UNIT=file_nr,FMT="(/,T2,A)") REPEAT("-",79)
      WRITE (UNIT=file_nr,FMT="(T2,A,T80,A)") "-","-"
      WRITE (UNIT=file_nr,FMT="(T2,A,T35,A,T80,A)") "-","TMC setting","-"
      WRITE (UNIT=file_nr,FMT="(T2,A,T80,A)") "-","-"
      WRITE (UNIT=file_nr,FMT="(T2,A)") REPEAT("-",79)

      WRITE (UNIT=file_nr,FMT="(T2,A,T35,A,T80,A)")"-","distribution of cores","-"
      WRITE(file_nr,FMT=fmt_my) plabel,"number of all working groups ", &
                      cp_to_string(tmc_env%tmc_comp_set%para_env_m_w%num_pe-1)
      WRITE(file_nr,FMT=fmt_my)plabel,"number of groups (ener|cc)",&
                      cp_to_string(tmc_env%tmc_comp_set%group_ener_nr)//" | "//&
                      cp_to_string(tmc_env%tmc_comp_set%group_cc_nr)
      WRITE(file_nr,FMT=fmt_my)plabel,"cores per group (ener|cc) ",&
                      cp_to_string(tmc_env%tmc_comp_set%group_ener_size)//" | "//&
                      cp_to_string(tmc_env%tmc_comp_set%group_cc_size)
      IF(ASSOCIATED(tmc_env%tmc_comp_set%para_env_m_ana)) &
         WRITE(file_nr,FMT=fmt_my)plabel,"Analysis groups ", &
           cp_to_string(tmc_env%tmc_comp_set%para_env_m_ana%num_pe-1)
      IF(SIZE(tmc_env%params%Temp(:)).LE.7) THEN
        WRITE(fmt_tmp,*)'(T2,A,"| ",A,T25,A56)'
        c_tmp = ""
        DO i = 1, SIZE(tmc_env%params%Temp(:))
          WRITE(c_tmp,FMT="(A,1000F8.2)")TRIM(c_tmp), tmc_env%params%Temp(i)
        END DO
        WRITE(file_nr,FMT=fmt_tmp)plabel,"Temperature(s) [K]", TRIM(c_tmp)
      ELSE
        WRITE(file_nr,FMT='(A,1000F8.2)')" "//plabel//"| Temperature(s) [K]",&
                        tmc_env%params%Temp(:)
      END IF
      WRITE(file_nr,FMT=fmt_my)plabel,"# of Monte Carlo Chain elements: ", &
                      cp_to_string(tmc_env%m_env%num_MC_elem)
      WRITE(file_nr,FMT=fmt_my)plabel,"exact potential input file:", &
                      TRIM(tmc_env%params%energy_inp_file)
      IF(tmc_env%params%NMC_inp_file.NE."") &
        WRITE(file_nr,FMT=fmt_my)plabel,"approximate potential input file:", &
                      TRIM(tmc_env%params%NMC_inp_file)
      IF(ANY(tmc_env%params%sub_box_size.GT.0.0_dp)) THEN
        WRITE(fmt_tmp,*)'(T2,A,"| ",A,T25,A56)'
        c_tmp = ""
        DO i = 1, SIZE(tmc_env%params%sub_box_size)
          WRITE(c_tmp,FMT="(A,1000F8.2)")TRIM(c_tmp), tmc_env%params%sub_box_size(i)*au2a
        END DO
        WRITE(file_nr,FMT=fmt_tmp)plabel,"Sub box size [A]", TRIM(c_tmp)
      END IF
      IF(tmc_env%params%pressure.GT.0.0_dp) &
        WRITE(file_nr,FMT=fmt_my)plabel,"Pressure [bar]: ", &
                      cp_to_string(tmc_env%params%pressure*au2bar)
      WRITE(file_nr,FMT=fmt_my)plabel,"Numbers of atoms/molecules moved "
      WRITE(file_nr,FMT=fmt_my)plabel,"  within one conf. change", &
                      cp_to_string(tmc_env%params%nr_elem_mv)
      WRITE (UNIT=file_nr,FMT="(/,T2,A)") REPEAT("-",79)
    END IF

  END SUBROUTINE tmc_print_params

END MODULE
