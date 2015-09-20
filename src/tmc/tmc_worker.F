!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief module contains the worker routine handling the communication and
!>        the calculation / creation of the configurations
!>        - WORKER these are all TMC cores, instead of master core 
!>          and maybe some idle cores
!>        - divided in groups, in every group exists group master
!>          - there can be two kind of groups, one for exact energy calculation
!>            and one calculating configurational change using an approximate 
!>            potential
!>        - Algorithm:
!>          - group master receive messages and decide what to do,
!>          - (if nessesary) broadcast of working task 
!>            to all other group members (needed for parallel CP2K)
!>          - process task, calculations of energy or configurational change
!>          - result, exist on group master, sent to master core
!>        Communication structure (master->worker, worker->master): 
!>        - message structure is defined in TMC message module
!> \par History
!>      11.2012 created [Mandes Schönherr]
!> \author Mandes
! *****************************************************************************

MODULE tmc_worker
  USE cell_types,                      ONLY: cell_copy,&
                                             cell_type,&
                                             init_cell
  USE cp_external_control,             ONLY: set_external_comm
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_result_methods,               ONLY: cp_results_erase,&
                                             put_results
  USE cp_result_types,                 ONLY: cp_result_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE f77_interface,                   ONLY: f_env_get_from_id,&
                                             f_env_type,&
                                             get_natom,&
                                             get_pos,&
                                             get_result_r1
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_get_natom
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE parallel_rng_types,              ONLY: get_rng_stream,&
                                             next_random_number,&
                                             set_rng_stream
  USE particle_list_types,             ONLY: particle_list_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE tmc_analysis,                    ONLY: analysis_init,&
                                             analysis_restart_print,&
                                             analysis_restart_read,&
                                             analyze_file_configurations,&
                                             do_tmc_analysis,&
                                             finalize_tmc_analysis
  USE tmc_analysis_types,              ONLY: tmc_ana_list_type
  USE tmc_calculations,                ONLY: calc_potential_energy
  USE tmc_messages,                    ONLY: bcast_group,&
                                             check_if_group_master,&
                                             communicate_atom_types,&
                                             master_comm_id,&
                                             recv_msg,&
                                             send_msg,&
                                             stop_whole_group,&
                                             tmc_message
  USE tmc_move_handle,                 ONLY: clear_move_probs,&
                                             prob_update,&
                                             select_random_move_type
  USE tmc_move_types,                  ONLY: mv_type_MD,&
                                             mv_type_NMC_moves
  USE tmc_moves,                       ONLY: change_pos
  USE tmc_stati,                       ONLY: &
       TMC_CANCELING_MESSAGE, TMC_CANCELING_RECEIPT, TMC_STATUS_CALCULATING, &
       TMC_STATUS_FAILED, TMC_STATUS_STOP_RECEIPT, &
       TMC_STATUS_WAIT_FOR_NEW_TASK, TMC_STATUS_WORKER_INIT, &
       TMC_STAT_ANALYSIS_REQUEST, TMC_STAT_ANALYSIS_RESULT, &
       TMC_STAT_APPROX_ENERGY_REQUEST, TMC_STAT_APPROX_ENERGY_RESULT, &
       TMC_STAT_ENERGY_REQUEST, TMC_STAT_ENERGY_RESULT, &
       TMC_STAT_INIT_ANALYSIS, TMC_STAT_MD_REQUEST, TMC_STAT_MD_RESULT, &
       TMC_STAT_NMC_REQUEST, TMC_STAT_NMC_RESULT, &
       TMC_STAT_SCF_STEP_ENER_RECEIVE, TMC_STAT_START_CONF_REQUEST, &
       TMC_STAT_START_CONF_RESULT, task_type_MC, task_type_ideal_gas
  USE tmc_tree_acceptance,             ONLY: acceptance_check
  USE tmc_tree_build,                  ONLY: allocate_new_sub_tree_node,&
                                             deallocate_sub_tree_node
  USE tmc_tree_types,                  ONLY: tree_type
  USE tmc_types,                       ONLY: allocate_tmc_atom_type,&
                                             tmc_atom_type,&
                                             tmc_env_type,&
                                             tmc_param_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tmc_worker'

  PUBLIC :: do_tmc_worker
  PUBLIC :: get_initial_conf, get_atom_kinds_and_cell

  INTEGER, PARAMETER :: DEBUG = 0

CONTAINS

! *****************************************************************************
!> \brief worker get tasks form master and fulfill them  
!> \param tmc_env structure for storing all the tmc parameters
!> \param ana_list ...
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE do_tmc_worker(tmc_env, ana_list)
    TYPE(tmc_env_type), POINTER              :: tmc_env
    TYPE(tmc_ana_list_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: ana_list

    CHARACTER(LEN=*), PARAMETER :: routineN = 'do_tmc_worker', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: c_tmp
    INTEGER                                  :: calc_stat, handle, i1, i2, &
                                                ierr, itmp, num_dim, work_stat
    INTEGER, DIMENSION(:), POINTER           :: ana_restart_conf
    LOGICAL                                  :: failure, flag, master
    TYPE(cp_para_env_type), POINTER          :: para_env_m_w
    TYPE(tree_type), POINTER                 :: conf

    failure = .FALSE.
    master  = .FALSE.
    i1 = -1
    i2 = -1
    NULLIFY(conf, para_env_m_w, ana_restart_conf)

    CPASSERT(ASSOCIATED(tmc_env))

    ! start the timing
    CALL timeset(routineN,handle)

    ! initialize
    IF(tmc_env%tmc_comp_set%group_nr.GT.0) THEN
      CPASSERT(ASSOCIATED(tmc_env%tmc_comp_set%para_env_sub_group))
      IF(tmc_env%w_env%env_id_ener.GT.0) THEN
        itmp = tmc_env%w_env%env_id_ener
      ELSE
        itmp = tmc_env%w_env%env_id_approx
      END IF
     
      CALL get_atom_kinds_and_cell(env_id=itmp, &
             atoms=tmc_env%params%atoms, cell=tmc_env%params%cell)
      para_env_m_w => tmc_env%tmc_comp_set%para_env_m_w
      master = check_if_group_master(tmc_env%tmc_comp_set%para_env_sub_group)
!      CALL get_mimal_distances()
    ELSE
      ! analysis group
      CPASSERT(ASSOCIATED(tmc_env%tmc_comp_set%para_env_m_ana))
      para_env_m_w => tmc_env%tmc_comp_set%para_env_m_ana
      master = .TRUE.
    END IF

      !-- GROUP MASTER only --------------
      ! get messages from master and handle them
      IF(master) THEN
        ! NOT the analysis group
        IF(tmc_env%tmc_comp_set%group_nr.GT.0) THEN
          IF(tmc_env%w_env%env_id_ener.GT.0) THEN
            itmp = tmc_env%w_env%env_id_ener
          ELSE
            itmp = tmc_env%w_env%env_id_approx
          END IF
          ! set the communicator in the external control for receiving exit tags 
          !  and sending additional information (e.g. the intermediate scf energies)
          IF(tmc_env%params%use_scf_energy_info) &
            CALL set_intermediate_info_comm(env_id=itmp, &
                  comm=tmc_env%tmc_comp_set%para_env_m_w%group)
          IF(tmc_env%params%SPECULATIVE_CANCELING) &
            CALL set_external_comm(comm=tmc_env%tmc_comp_set%para_env_m_w%group, &
                   in_external_master_id=MASTER_COMM_ID, &
                   in_exit_tag=TMC_CANCELING_MESSAGE)
        END IF
        !-- WORKING LOOP --!
        master_work_time: DO
          work_stat = TMC_STATUS_WAIT_FOR_NEW_TASK
          ! -- receive message from master
          ! check for new task (wait for it)
          itmp = MASTER_COMM_ID
          CALL tmc_message(msg_type=work_stat, send_recv=recv_msg, &
                           dest=itmp, &
                           para_env=para_env_m_w, &
                           result_count=ana_restart_conf,&
                           tmc_params=tmc_env%params, elem=conf)

          IF(DEBUG.GE.1.AND.work_stat.NE.TMC_STATUS_WAIT_FOR_NEW_TASK) &
            WRITE(tmc_env%w_env%io_unit,*)"worker: group master of group ", &
              tmc_env%tmc_comp_set%group_nr,"got task ", work_stat
          calc_stat = TMC_STATUS_CALCULATING
          SELECT CASE(work_stat)
          CASE(TMC_STATUS_WAIT_FOR_NEW_TASK)
          CASE(TMC_STATUS_WORKER_INIT)
            CALL init_cell(cell=tmc_env%params%cell)
            itmp = bcast_group
            CALL tmc_message(msg_type=work_stat, send_recv=send_msg, &
                             dest=itmp, &
                             para_env=tmc_env%tmc_comp_set%para_env_sub_group,&
                             tmc_params=tmc_env%params)
          CASE(TMC_CANCELING_MESSAGE)
            work_stat = TMC_CANCELING_RECEIPT
            itmp = MASTER_COMM_ID
            CALL tmc_message(msg_type=work_stat, send_recv=send_msg, &
                             dest=itmp, &
                           para_env=para_env_m_w, &
                           tmc_params=tmc_env%params)
          CASE(TMC_STATUS_FAILED)
            IF(DEBUG.GE.1) &
              WRITE(tmc_env%w_env%io_unit,*) "master worker of group", &
                tmc_env%tmc_comp_set%group_nr, " exit work time."
            EXIT master_work_time
          !-- group master read the CP2K input file, and write data to master
          CASE(TMC_STAT_START_CONF_REQUEST)
            IF(tmc_env%w_env%env_id_ener.GT.0) THEN
              itmp = tmc_env%w_env%env_id_ener
            ELSE
              itmp = tmc_env%w_env%env_id_approx
            END IF
            CALL get_initial_conf(tmc_params=tmc_env%params,init_conf=conf,&
                                  env_id=itmp)
            ! send start configuration back to master
            work_stat = TMC_STAT_START_CONF_RESULT
            itmp = MASTER_COMM_ID
            CALL tmc_message(msg_type=work_stat, send_recv=send_msg, &
                             dest=itmp, &
                             para_env=para_env_m_w, &
                             tmc_params=tmc_env%params, elem=conf, &
                             wait_for_message=.TRUE.)
            
            IF(ASSOCIATED(tmc_env%tmc_comp_set%para_env_m_first_w))&
              CALL communicate_atom_types(atoms=tmc_env%params%atoms, &
                     source=1, &
                     para_env=tmc_env%tmc_comp_set%para_env_m_first_w)
          !-- calculate the approximate energy
          CASE(TMC_STAT_APPROX_ENERGY_REQUEST)
            CPASSERT(tmc_env%w_env%env_id_approx.GT.0)
            itmp = bcast_group
            !-- DISTRIBUTING WORK (group master) to all other group members
            CALL tmc_message(msg_type=work_stat, send_recv=send_msg, &
                             dest=itmp, &
                             para_env=tmc_env%tmc_comp_set%para_env_sub_group,&
                             tmc_params=tmc_env%params, elem=conf)
            CALL calc_potential_energy(conf=conf, &
                                       env_id=tmc_env%w_env%env_id_approx, &
                                       exact_approx_pot=.FALSE., &
                                       tmc_env=tmc_env) 
            work_stat = TMC_STAT_APPROX_ENERGY_RESULT
            itmp = MASTER_COMM_ID
            CALL tmc_message(msg_type=work_stat, send_recv=send_msg, &
                             dest=itmp, &
                             para_env=para_env_m_w, &
                             tmc_params=tmc_env%params, elem=conf)
          ! -- Nested Monte Carlo routines
          CASE(TMC_STAT_MD_REQUEST, TMC_STAT_NMC_REQUEST)
            CALL clear_move_probs(tmc_env%params%nmc_move_types)
            itmp = bcast_group
            CALL tmc_message(msg_type=work_stat, send_recv=send_msg, &
                             dest=itmp, &
                             para_env=tmc_env%tmc_comp_set%para_env_sub_group,&
                             tmc_params=tmc_env%params, elem=conf)
            !-- collective calculation for MD/NMC steps
            IF(work_stat.EQ.TMC_STAT_NMC_REQUEST)THEN
              !-- calculate MD steps, in case of 2 different potentials do nested Monte Carlo
              CALL nested_markov_chain_MC(conf=conf, &
                     env_id=tmc_env%w_env%env_id_approx,&
                     tmc_env=tmc_env, calc_status=calc_stat)
            ELSEIF(work_stat.EQ.TMC_STAT_MD_REQUEST) THEN
              !TODO Hybrid MC routine
              CALL cp_assert(.FALSE.,&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "there is no Hybrid MC implemented yet.",&
                     only_ionode=.TRUE.)

            ELSE
              CALL cp_assert(.FALSE.,&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "unknown task type for workers.",&
                     only_ionode=.TRUE.)
            END IF
            !-- in case of cancelation send receipt
            itmp = MASTER_COMM_ID
            CALL tmc_message(msg_type=calc_stat, send_recv=recv_msg, &
                             dest=itmp, &
                             para_env=para_env_m_w, &
                             tmc_params=tmc_env%params, &
                             success=flag)
            SELECT CASE(calc_stat)
            CASE(TMC_STATUS_CALCULATING)
              SELECT CASE(work_stat)
              CASE(TMC_STAT_MD_REQUEST)
                work_stat = TMC_STAT_MD_RESULT
              CASE(TMC_STAT_NMC_REQUEST)
                work_stat = TMC_STAT_NMC_RESULT
              CASE DEFAULT
                CALL cp_assert(.FALSE., cp_failure_level,cp_assertion_failed,&
                       routineP,&
                       "unknown work status after possible NMC subgroup "//&
                       "cancelation, work_stat="//cp_to_string(work_stat))
              END SELECT
            CASE(TMC_CANCELING_MESSAGE)
              work_stat = TMC_CANCELING_RECEIPT
            CASE DEFAULT
              CALL cp_assert(.FALSE., cp_failure_level,cp_assertion_failed,&
                       routineP,&
                       "unknown calc status before sending NMC result "//&
                       cp_to_string(calc_stat))
            END SELECT
            ! send message back to master
            itmp = MASTER_COMM_ID
            CALL tmc_message(msg_type=work_stat, send_recv=send_msg, &
                             dest=itmp, &
                             para_env=para_env_m_w, &
                             tmc_params=tmc_env%params, elem=conf)
          CASE(TMC_STAT_ENERGY_REQUEST)
            CPASSERT(tmc_env%w_env%env_id_ener.GT.0)
            !-- DISTRIBUTING WORK (group master) to all other group members
            itmp = bcast_group
            CALL tmc_message(msg_type=work_stat, send_recv=send_msg, &
                             dest=itmp, &
                             para_env=tmc_env%tmc_comp_set%para_env_sub_group,&
                             tmc_params=tmc_env%params, elem=conf)

            CALL calc_potential_energy(conf=conf, &
                                       env_id=tmc_env%w_env%env_id_ener, &
                                       exact_approx_pot=.TRUE., &
                                       tmc_env=tmc_env)
            !-- in case of cancelation send receipt
            flag = .FALSE.
            itmp = MASTER_COMM_ID
            CALL tmc_message(msg_type=calc_stat, send_recv=recv_msg, &
                             dest=itmp, &
                             para_env=para_env_m_w,&
                             tmc_params=tmc_env%params, success=flag)
            SELECT CASE(calc_stat)
            CASE(TMC_STATUS_CALCULATING)
              SELECT CASE(work_stat)
              CASE(TMC_STAT_ENERGY_REQUEST)
                work_stat = TMC_STAT_ENERGY_RESULT
                !-- if nessesary get the exact dipoles (for e.g. quantum potential)
                IF(tmc_env%params%print_dipole) THEN
                  c_tmp = "[DIPOLE]"
                  CALL get_result_r1(env_id=tmc_env%w_env%env_id_ener, &
                                     description=c_tmp,N=3, RESULT=conf%dipole,&
                                     res_exist=flag, ierr=ierr)
                  IF(.NOT.flag) tmc_env%params%print_dipole = .FALSE.
                  ! TODO maybe let run with the changed option, but inform user properly
                  CALL cp_assert(flag, cp_failure_level,cp_assertion_failed,&
                       routineP,&
                       "TMC: The requested dipoles are not porvided by the "//&
                       "force environment.")
                END IF
              CASE DEFAULT
                CALL cp_assert(.FALSE., cp_failure_level,cp_assertion_failed,&
                       routineP,&
                       "energy worker should handle unknown stat "//&
                       cp_to_string(work_stat))
              END SELECT
            CASE(TMC_CANCELING_MESSAGE)
               work_stat = TMC_CANCELING_RECEIPT
            CASE DEFAULT
               CALL cp_assert(.FALSE., cp_failure_level,cp_assertion_failed,&
                       routineP,&
                       "worker while energy calc is in unknown state "//&
                       cp_to_string(work_stat))
            END SELECT

            !-- send information back to master
            IF(DEBUG.GE.1) &
              WRITE(tmc_env%w_env%io_unit,*)"worker group ", &
                    tmc_env%tmc_comp_set%group_nr,&
                    "calculations done, send result energy", conf%potential
            itmp = MASTER_COMM_ID
            CALL tmc_message(msg_type=work_stat, send_recv=send_msg,&
                             dest=itmp, &
                             para_env=para_env_m_w,&
                             tmc_params=tmc_env%params, elem=conf)
          CASE(TMC_STAT_INIT_ANALYSIS)
            CPASSERT(ASSOCIATED(ana_restart_conf))
            CPASSERT(SIZE(ana_restart_conf).EQ.tmc_env%params%nr_temp)
            CPASSERT(PRESENT(ana_list))
            CPASSERT(ASSOCIATED(ana_list))
            itmp = MASTER_COMM_ID
            CALL communicate_atom_types(atoms=tmc_env%params%atoms, &
                   source=itmp, para_env=tmc_env%tmc_comp_set%para_env_m_ana)
            
            num_dim = SIZE(conf%pos)
            DO itmp=1, tmc_env%params%nr_temp
              ! do not forget to nullify the pointer at the end, deallcoated at tmc_env%params
              ana_list(itmp)%temp%temperature = tmc_env%params%Temp(itmp)
              ana_list(itmp)%temp%atoms       => tmc_env%params%atoms 
              ana_list(itmp)%temp%cell        => tmc_env%params%cell
!              ana_list(itmp)%temp%io_unit     = tmc_env%w_env%io_unit
              
              CALL analysis_init(ana_env=ana_list(itmp)%temp, nr_dim=num_dim) 
              ana_list(itmp)%temp%print_test_output = tmc_env%params%print_test_output
              IF(.NOT.ASSOCIATED(conf)) &
                CALL allocate_new_sub_tree_node(tmc_params=tmc_env%params, &
                                            next_el=conf, nr_dim=num_dim)
              CALL analysis_restart_read(ana_env=ana_list(itmp)%temp, &
                                         elem=conf)
              !check if we have the read the file
              flag = .FALSE.
              IF((.NOT.ASSOCIATED(ana_list(itmp)%temp%last_elem)) .AND. &
                 ana_restart_conf(itmp).GT.0) THEN
                flag = .TRUE.
                i1 = 0
                i2 = ana_restart_conf(itmp)
                CALL cp_assert(.FALSE., &
                       cp_warning_level,cp_assertion_failed,&
                       routineP,&
                       "analysis old trajectory up to "//&
                       "elem "//cp_to_string(ana_restart_conf(itmp))//&
                       ". Read trajectory file.")
              ELSE IF(ASSOCIATED(ana_list(itmp)%temp%last_elem))THEN
                IF(.NOT.(ana_list(itmp)%temp%last_elem%nr.EQ.ana_restart_conf(itmp))) THEN
                  flag = .TRUE.
                  i1 = ana_list(itmp)%temp%last_elem%nr
                  i2 = ana_restart_conf(itmp)
                  CALL cp_assert(.FALSE., &
                         cp_warning_level,cp_assertion_failed,&
                         routineP,&
                         "analysis restart with the incorrect configuration "//&
                         "TMC "//cp_to_string(ana_restart_conf(itmp))//&
                         " ana "//cp_to_string(ana_list(itmp)%temp%last_elem%nr)//&
                         ". REread trajectory file.")
                END IF
              END IF
              IF(flag) THEN
                CALL analyze_file_configurations(start_id=i1, &
                                                 end_id=i2, &
                                                 ana_env=ana_list(itmp)%temp, &
                                                 tmc_params=tmc_env%params)
              END IF
            END DO 
          CASE(TMC_STAT_ANALYSIS_REQUEST)
            CPASSERT(PRESENT(ana_list))
            CPASSERT(ASSOCIATED(ana_list(conf%sub_tree_nr)%temp))
            CALL do_tmc_analysis(elem=conf, &
                   ana_env=ana_list(conf%sub_tree_nr)%temp)
            work_stat = TMC_STAT_ANALYSIS_RESULT
            itmp = MASTER_COMM_ID
            CALL tmc_message(msg_type=work_stat, send_recv=send_msg,&
                             dest=itmp, &
                             para_env=para_env_m_w,&
                             tmc_params=tmc_env%params, elem=conf)
          CASE DEFAULT
            CALL cp_assert(.FALSE., cp_failure_level,cp_assertion_failed,&
                       routineP,&
                       "worker received unknown message task type "//&
                       cp_to_string(work_stat))
          END SELECT

          IF(DEBUG.GE.1.AND.work_stat.NE.TMC_STATUS_WAIT_FOR_NEW_TASK) &
            WRITE(tmc_env%w_env%io_unit,*)"worker: group ", &
              tmc_env%tmc_comp_set%group_nr, &
              "send back status:", work_stat
          IF(ASSOCIATED(conf)) &
            CALL deallocate_sub_tree_node(tree_elem=conf)
        END DO master_work_time
      !-- every other group paricipants----------------------------------------
      ELSE
        worker_work_time: DO
          work_stat = TMC_STATUS_WAIT_FOR_NEW_TASK
          flag = .FALSE.
          itmp = bcast_group
          CALL tmc_message(msg_type=work_stat, send_recv=recv_msg, &
                           dest=itmp, &
                           para_env=tmc_env%tmc_comp_set%para_env_sub_group,&
                           tmc_params=tmc_env%params, elem=conf)
          calc_stat = TMC_STATUS_CALCULATING
          SELECT CASE(work_stat)
          CASE(TMC_STATUS_WORKER_INIT)
            CALL init_cell(cell=tmc_env%params%cell)
          CASE(TMC_CANCELING_MESSAGE)
          ! error message
          CASE(TMC_STATUS_FAILED)
            EXIT worker_work_time
          ! all group members have to calculate the (MD potential) energy together
          CASE(TMC_STAT_START_CONF_RESULT)
            CPASSERT(tmc_env%w_env%env_id_approx.GT.0)
            !-- collective calculation of the potential energy of MD potential
            SELECT CASE(tmc_env%params%task_type)
            CASE(task_type_MC, task_type_ideal_gas)
              IF(tmc_env%params%NMC_inp_file.NE."") THEN
                conf%box_scale(:) = 1.0_dp
                CALL calc_potential_energy(conf=conf, &
                                           env_id=tmc_env%w_env%env_id_approx, &
                                           exact_approx_pot=.FALSE., &
                                           tmc_env=tmc_env)
              END IF
            CASE DEFAULT
              CALL cp_assert(.FALSE., cp_failure_level,cp_assertion_failed,&
                             routineP,&
                             "unknown task_type for participants in "//&
                             "START_CONF_RESULT request ")
            END SELECT
          !-- HMC - calculating MD steps
          CASE(TMC_STAT_NMC_REQUEST, TMC_STAT_MD_REQUEST)
            !-- collective calculation for MD/NMC steps
            IF(work_stat.EQ.TMC_STAT_NMC_REQUEST)THEN
              !-- calculate MD steps, in case of 2 different potentials do nested Monte Carlo
              CALL nested_markov_chain_MC(conf=conf, &
                     env_id=tmc_env%w_env%env_id_approx,&
                     tmc_env=tmc_env, calc_status=calc_stat)
            ELSEIF(work_stat.EQ.TMC_STAT_MD_REQUEST) THEN
              !TODO Hybrid MC routine
              CALL cp_assert(.FALSE.,&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "there is no Hybrid MC implemented yet.",&
                     only_ionode=.TRUE.)

            ELSE
              CALL cp_assert(.FALSE.,&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "unknown task type for workers.",&
                     only_ionode=.TRUE.)
            END IF
          !-- energy calculations
          CASE(TMC_STAT_APPROX_ENERGY_REQUEST)
            !--- do calculate energy
            CPASSERT(tmc_env%w_env%env_id_approx.GT.0)
            CALL calc_potential_energy(conf=conf, &
                                       env_id=tmc_env%w_env%env_id_approx, &
                                       exact_approx_pot=.FALSE., &
                                       tmc_env=tmc_env)
          CASE(TMC_STAT_ENERGY_REQUEST)
            !--- do calculate energy
            CPASSERT(tmc_env%w_env%env_id_ener.GT.0)
            CALL calc_potential_energy(conf=conf, &
                                       env_id=tmc_env%w_env%env_id_ener, &
                                       exact_approx_pot=.TRUE., &
                                       tmc_env=tmc_env)
          CASE DEFAULT
            CALL cp_assert(.FALSE., cp_failure_level,cp_assertion_failed,&
                       routineP,&
                       "group participant got unknown working type "//&
                       cp_to_string(work_stat))
          END SELECT
          IF(ASSOCIATED(conf)) &
            CALL deallocate_sub_tree_node(tree_elem=conf)
        END DO worker_work_time
      END IF
      ! --------------------------------------------------------------------
      ! finalizing analysis, writing files etc.
      IF(ASSOCIATED(tmc_env%tmc_comp_set%para_env_m_ana))THEN
        DO itmp=1, tmc_env%params%nr_temp
          CALL analysis_restart_print(ana_env=ana_list(itmp)%temp)
          IF(ASSOCIATED(conf)) &
            CALL deallocate_sub_tree_node(tree_elem=ana_list(itmp)%temp%last_elem)
          CALL finalize_tmc_analysis(ana_list(itmp)%temp)
        END DO
      END IF
      !-- stopping and finalizing
      ! sending back receipt for stopping
      IF(master) THEN
        ! NOT the analysis group
        IF(tmc_env%tmc_comp_set%group_nr.GT.0) THEN
          ! remove the communicator in the external control for receiving exit tags 
          !  and sending additional information (e.g. the intermediate scf energies)
          IF(tmc_env%params%use_scf_energy_info) THEN
            IF(tmc_env%w_env%env_id_ener.GT.0) THEN
              itmp = tmc_env%w_env%env_id_ener
            ELSE
              itmp = tmc_env%w_env%env_id_approx
            END IF
            CALL remove_intermediate_info_comm(env_id=itmp)
          END IF
        END IF
        IF(ASSOCIATED(tmc_env%tmc_comp_set%para_env_sub_group)) &
          CALL stop_whole_group(para_env=tmc_env%tmc_comp_set%para_env_sub_group,&
                                tmc_params=tmc_env%params)
        
          work_stat = TMC_STATUS_STOP_RECEIPT
          itmp = MASTER_COMM_ID
          CALL tmc_message(msg_type=work_stat, send_recv=send_msg, dest=itmp,&
                           para_env=para_env_m_w,&
                           tmc_params=tmc_env%params)
      ELSE IF(ASSOCIATED(tmc_env%tmc_comp_set%para_env_sub_group))THEN
        work_stat = TMC_STATUS_STOP_RECEIPT
        itmp = MASTER_COMM_ID
        CALL tmc_message(msg_type=work_stat, send_recv=send_msg, dest=itmp,&
                         para_env=tmc_env%tmc_comp_set%para_env_sub_group,&
                         tmc_params=tmc_env%params)
      END IF

      IF(DEBUG.GE.5)&
        WRITE(tmc_env%w_env%io_unit,*)"worker ", &
          tmc_env%tmc_comp_set%para_env_sub_group%mepos, "of group ",&
          tmc_env%tmc_comp_set%group_nr,"stops working!!!!!!!!!!!!!!!!!!"

    IF(PRESENT(ana_list)) THEN
      DO itmp=1, tmc_env%params%nr_temp
        ana_list(itmp)%temp%atoms => NULL() 
        ana_list(itmp)%temp%cell  => NULL()
      END DO
    END IF
    IF(ASSOCIATED(conf)) &
      CALL deallocate_sub_tree_node(tree_elem=conf)
    IF(ASSOCIATED(ana_restart_conf)) DEALLOCATE(ana_restart_conf)

    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE do_tmc_worker

! *****************************************************************************
!> \brief Nested Monte Carlo (NMC), do several Markov Chain Monte Carlo steps 
!>        usually using the approximate potential, could be also Hybrid MC.
!>        The amount of steps are predefinded by the user, but should be huge 
!>        enough to reach the equilibrium state for this potential
!> \param conf ...
!> \param env_id ...
!> \param tmc_env ...
!> \param calc_status ...
!> \param 
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE nested_markov_chain_MC(conf, env_id, tmc_env, calc_status)
    TYPE(tree_type), POINTER                 :: conf
    INTEGER, INTENT(IN)                      :: env_id
    TYPE(tmc_env_type), POINTER              :: tmc_env
    INTEGER, INTENT(OUT)                     :: calc_status

    CHARACTER(LEN=*), PARAMETER :: routineN = 'nested_markov_chain_MC', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: comm_dest, handle, substeps
    LOGICAL                                  :: accept, change_rejected, flag
    REAL(KIND=dp)                            :: rnd_nr
    TYPE(tree_type), POINTER                 :: last_acc_conf

    NULLIFY(last_acc_conf)
    
    CPASSERT(ASSOCIATED(tmc_env))
    CPASSERT(ASSOCIATED(tmc_env%params))
    CPASSERT(ASSOCIATED(tmc_env%tmc_comp_set))
    CPASSERT(ASSOCIATED(tmc_env%rng_stream))
    CPASSERT(ASSOCIATED(conf))
    CPASSERT(conf%temp_created.GT.0)
    CPASSERT(conf%temp_created.LE.tmc_env%params%nr_temp)
    CPASSERT(env_id.GT.0)

    ! start the timing
    CALL timeset(routineN,handle)

    CALL allocate_new_sub_tree_node(tmc_params=tmc_env%params, &
           next_el=last_acc_conf, nr_dim=SIZE(conf%pos))

    last_acc_conf%pos       = conf%pos
    last_acc_conf%box_scale = conf%box_scale

    ! energy of the last accepted configuration
    CALL calc_potential_energy(conf=last_acc_conf, &
           env_id=tmc_env%w_env%env_id_approx,exact_approx_pot=.FALSE.,&
           tmc_env=tmc_env)

    NMC_steps: DO substeps=1, INT(tmc_env%params%move_types%mv_size(mv_type_NMC_moves,1))
      ! check for canceling message
      IF(ASSOCIATED(tmc_env%tmc_comp_set%para_env_m_w))THEN
        flag = .FALSE.
        comm_dest = MASTER_COMM_ID
        ! check for new canceling message
        CALL tmc_message(msg_type=calc_status, send_recv=recv_msg, &
                         dest=comm_dest, &
                         para_env=tmc_env%tmc_comp_set%para_env_m_w, &
                         tmc_params=tmc_env%params, success=flag)
      END IF
      comm_dest = bcast_group
      CALL tmc_message(msg_type=calc_status, send_recv=send_msg, &
                       dest=comm_dest, &
                       para_env=tmc_env%tmc_comp_set%para_env_sub_group, &
                       tmc_params=tmc_env%params)
      SELECT CASE(calc_status)
      CASE(TMC_STATUS_CALCULATING)
        ! keep on working
      CASE(TMC_CANCELING_MESSAGE)
        ! nothing to do, because calculation CANCELING, exit with cancel status
        EXIT NMC_steps
      CASE DEFAULT
        CALL cp_assert(.FALSE., cp_failure_level,cp_assertion_failed,routineP,&
                         "unknown status "//cp_to_string(calc_status)//&
                         "in the NMC routine, expect only caneling status. ")
      END SELECT

      ! set move type
      CALL set_rng_stream(rng_stream=tmc_env%rng_stream,&
             bg=conf%rng_seed(:,:,1), cg=conf%rng_seed(:,:,2),&
             ig=conf%rng_seed(:,:,3))
      conf%move_type = select_random_move_type(&
                         move_types=tmc_env%params%nmc_move_types, &
                         rnd=next_random_number(tmc_env%rng_stream))
      CALL get_rng_stream(rng_stream=tmc_env%rng_stream,&
             bg=conf%rng_seed(:,:,1), cg=conf%rng_seed(:,:,2),&
             ig=conf%rng_seed(:,:,3))

      ! do move
      CALL change_pos(tmc_params=tmc_env%params, &
                      move_types=tmc_env%params%nmc_move_types, &
                      rng_stream=tmc_env%rng_stream,&
                      elem=conf, mv_conf=1, new_subbox=.FALSE., &
                      move_rejected=change_rejected)
      ! for Hybrid MC the change_pos is only velocity change, 
      !   the actual MD step hast to be done in this module for communication reason
      IF(conf%move_type .EQ. mv_type_MD)THEN
        !TODO implement the MD part
        !CALL calc_MD_step(...)
        !CALL calc_calc_e_kin(...)
        CALL cp_assert(.FALSE., cp_failure_level,cp_assertion_failed,routineP,&
                       "Hybrid MC is not implemented yet, "//&
                       "(no MD section in TMC yet). ")
      END IF

      ! update the subbox acceptance probabilities
      CALL prob_update(move_types=tmc_env%params%nmc_move_types, elem=conf, &
                       acc=.NOT.change_rejected, subbox=.TRUE., &
                       prob_opt=tmc_env%params%esimate_acc_prob)

      ! calculate potential energy if neccessary
      IF(.NOT.change_rejected) THEN
        CALL calc_potential_energy(conf=conf, &
               env_id=tmc_env%w_env%env_id_approx, exact_approx_pot=.FALSE., &
               tmc_env=tmc_env)
      ELSE
         conf%e_pot_approx = HUGE(conf%e_pot_approx)
      END IF

      !check NMC step
      CALL set_rng_stream(rng_stream=tmc_env%rng_stream, &
             bg=conf%rng_seed(:,:,1), cg=conf%rng_seed(:,:,2), &
             ig=conf%rng_seed(:,:,3))
      rnd_nr = next_random_number(tmc_env%rng_stream)
      CALL get_rng_stream(rng_stream=tmc_env%rng_stream, &
             bg=conf%rng_seed(:,:,1), cg=conf%rng_seed(:,:,2), &
             ig=conf%rng_seed(:,:,3))

      IF(.NOT.change_rejected) THEN
        CALL acceptance_check(tree_element=conf, parent_element=last_acc_conf,&
               tmc_params=tmc_env%params, &
               temperature=tmc_env%params%Temp(conf%temp_created), &
               diff_pot_check=.FALSE.,&
               accept=accept, approx_ener=.TRUE., rnd_nr=rnd_nr)
      ELSE
        accept=.FALSE.
      END IF
      ! update the NMC accpetance per move
      CALL prob_update(move_types=tmc_env%params%nmc_move_types, elem=conf, &
             acc=accept, prob_opt=tmc_env%params%esimate_acc_prob)

      ! update last accepted configuration or actual configuration
      IF(accept.AND.(.NOT.change_rejected))THEN
        last_acc_conf%pos            = conf%pos
        last_acc_conf%vel            = conf%vel
        last_acc_conf%e_pot_approx   = conf%e_pot_approx
        last_acc_conf%ekin           = conf%ekin
        last_acc_conf%ekin_before_md = conf%ekin_before_md
        last_acc_conf%box_scale      = conf%box_scale
      ELSE
        conf%pos       = last_acc_conf%pos
        conf%vel       = last_acc_conf%vel
        conf%box_scale = last_acc_conf%box_scale
      END IF
    END DO NMC_steps

    ! result values of Nested Monte Carlo (NMC) steps 
    !   regard that the calculated potential energy is the one of the approximated potential
    conf%pos            = last_acc_conf%pos
    conf%vel            = last_acc_conf%vel
    conf%e_pot_approx   = last_acc_conf%e_pot_approx
    conf%potential      = 0.0_dp
    conf%ekin           = last_acc_conf%ekin
    conf%ekin_before_md = last_acc_conf%ekin_before_md

    CALL deallocate_sub_tree_node(tree_elem=last_acc_conf)

    ! end the timing 
    CALL timestop(handle)
  END SUBROUTINE nested_markov_chain_MC

! *****************************************************************************
!> \brief get the initial confuguration (pos,...) 
!> \param tmc_params ...
!> \param init_conf the structure the data should be stored
!> force_env 
!> \param env_id ...
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE get_initial_conf(tmc_params, init_conf, env_id)
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(tree_type), POINTER                 :: init_conf
    INTEGER                                  :: env_id

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_initial_conf', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, mol, ndim, &
                                                nr_atoms
    LOGICAL                                  :: failure
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(mol_new_list_type), POINTER         :: molecule_new

    failure = .FALSE.
    CPASSERT(.NOT.ASSOCIATED(init_conf))

    ! start the timing
    CALL timeset(routineN,handle)

    ! get positions
    CALL get_natom(env_id=env_id, n_atom=nr_atoms, ierr=ierr)
    CPASSERT(ierr.EQ.0)
    ndim  = 3*nr_atoms
    CALL allocate_new_sub_tree_node(tmc_params=tmc_params, &
                                    next_el=init_conf, nr_dim=ndim)
    CALL get_pos(env_id=env_id, pos=init_conf%pos, n_el=SIZE(init_conf%pos), &
                 ierr=ierr)

    ! get the molecule info 
    CALL f_env_get_from_id(env_id,f_env)
    CALL force_env_get(f_env%force_env,subsys=subsys)

    CALL cp_subsys_get(subsys=subsys, molecules_new=molecule_new)
    loop_mol: DO mol=1, SIZE(molecule_new%els(:))
      init_conf%mol(molecule_new%els(mol)%first_atom:&
                    molecule_new%els(mol)%last_atom) = mol
    END DO loop_mol

    ! end the timing
    CALL timestop(handle)

  END SUBROUTINE get_initial_conf

! *****************************************************************************
!> \brief get the pointer to the atoms, for easy handling
!> \param env_id ...
!> \param atoms pointer to atomic_kind 
!> \param cell ...
!> \author Mandes 01.2013
! *****************************************************************************
  SUBROUTINE get_atom_kinds_and_cell(env_id, atoms, cell)
    INTEGER                                  :: env_id
    TYPE(tmc_atom_type), DIMENSION(:), &
      POINTER                                :: atoms
    TYPE(cell_type), POINTER                 :: cell

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_atom_kinds_and_cell', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iparticle, nr_atoms, &
                                                nunits_tot
    LOGICAL                                  :: failure
    TYPE(cell_type), POINTER                 :: cell_tmp
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(particle_list_type), POINTER        :: particles

    failure = .FALSE.
    NULLIFY (f_env, subsys, particles)
    nr_atoms = 0

    CPASSERT(env_id.GT.0)
    CPASSERT(.NOT.ASSOCIATED(atoms))
    CPASSERT(.NOT.ASSOCIATED(cell))

    CALL f_env_get_from_id(env_id,f_env)
    nr_atoms = force_env_get_natom(f_env%force_env)
    CALL force_env_get(f_env%force_env,subsys=subsys,cell=cell_tmp)
    ALLOCATE(cell)
    CALL cell_copy(cell_in=cell_tmp, cell_out=cell)

    !get atom kinds
    CALL allocate_tmc_atom_type(atoms, nr_atoms)
    CALL cp_subsys_get(subsys, particles=particles)
    nunits_tot=SIZE(particles%els(:))
    IF(nunits_tot .GT. 0) THEN
       DO iparticle=1, nunits_tot
         atoms(iparticle)%name = particles%els(iparticle)%atomic_kind%name
         atoms(iparticle)%mass = particles%els(iparticle)%atomic_kind%mass
       END DO
       CPASSERT(iparticle-1.EQ.nr_atoms)
    ENDIF
  END SUBROUTINE get_atom_kinds_and_cell

! *****************************************************************************
!> \brief set the communicator in the SCF environment 
!>        to receive the intermediate energies on the (global) master side
!> \param comm the master-worker communicator
!> \param env_id the ID of the related force environment
!> \author Mandes 10.2013
! *****************************************************************************
  SUBROUTINE set_intermediate_info_comm(comm, env_id)
    INTEGER, INTENT(IN)                      :: comm
    INTEGER                                  :: env_id

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_intermediate_info_comm', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: description
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(3)              :: values
    TYPE(cp_result_type), POINTER            :: results
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(f_env_type), POINTER                :: f_env

    failure = .FALSE.
    NULLIFY(results,subsys)
    CPASSERT(env_id.GT.0)

    CALL f_env_get_from_id(env_id,f_env)

    CPASSERT(ASSOCIATED(f_env))
    CPASSERT(ASSOCIATED(f_env%force_env))
    CALL cp_assert(ASSOCIATED(f_env%force_env%qs_env), &
                   cp_failure_level,cp_assertion_failed,routineP,&
                   "the intermediate SCF energy request can not be set "//&
                   "employing this force environment! ") 

    ! set the information
    values(1) = REAL(comm, KIND=dp)
    values(2) = REAL(MASTER_COMM_ID, KIND=dp)
    values(3) = REAL(TMC_STAT_SCF_STEP_ENER_RECEIVE, KIND=dp)
    description = "[EXT_SCF_ENER_COMM]"

    ! set the communicator information in the qs_env result container
    CALL force_env_get(f_env%force_env, subsys=subsys)
    CALL cp_subsys_get(subsys, results=results)
    CALL put_results(results, description=description, values=values)
  END SUBROUTINE set_intermediate_info_comm

! *****************************************************************************
!> \brief set the communicator in the SCF environment 
!>        to receive the intermediate energies on the (global) master side
!> \param env_id the ID of the related force environment
!> \author Mandes 10.2013
! *****************************************************************************
  SUBROUTINE remove_intermediate_info_comm(env_id)
    INTEGER                                  :: env_id

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'remove_intermediate_info_comm', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: description
    LOGICAL                                  :: failure
    TYPE(cp_result_type), POINTER            :: results
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(f_env_type), POINTER                :: f_env

    failure = .FALSE.
    NULLIFY(subsys,results)
    CPASSERT(env_id.GT.0)

    CALL f_env_get_from_id(env_id,f_env)

    CPASSERT(ASSOCIATED(f_env))
    CPASSERT(ASSOCIATED(f_env%force_env))
    CALL cp_assert(ASSOCIATED(f_env%force_env%qs_env), &
                   cp_failure_level,cp_assertion_failed,routineP,&
                   "the SCF intermediate energy communicator can not be "//&
                   "removed! ")

    description = "[EXT_SCF_ENER_COMM]"

    ! set the communicator information in the qs_env result container
    CALL force_env_get(f_env%force_env, subsys=subsys)
    CALL cp_subsys_get(subsys, results=results)
    CALL cp_results_erase(results, description=description)
  END SUBROUTINE remove_intermediate_info_comm


!! *****************************************************************************
!!> \brief get the pointer to the minimal distances 
!!>        to check for overlapping atoms
!!> \param rlist_lowsq minimal distanced depending on atom types
!!> \param error variable to control error logging, stopping,...
!!>        see module cp_error_handling
!!> \author Mandes 03.2013
!! *****************************************************************************
!  SUBROUTINE get_mimal_distances(env_id, rlist_lowsq)
!    INTEGER                                  :: env_id
!    REAL(KIND=dp), DIMENSION(:, :), &
!      OPTIONAL, POINTER                      :: rlist_lowsq
!
!    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_mimal_distances', &
!      routineP = moduleN//':'//routineN
!
!    LOGICAL                                  :: failure
!    TYPE(f_env_type), POINTER                :: f_env
!
!    CPPrecondition(env_id.GT.0,cp_failure_level,routineP,failure)
!    CPPrecondition(.NOT.ASSOCIATED(rlist_lowsq),cp_failure_level,routineP,failure)
!
!    CALL f_env_get_from_id(env_id,f_env)
!    IF (.NOT.failure) THEN
!      IF(ASSOCIATED(f_env%force_env)) THEN
!        IF(ASSOCIATED(f_env%force_env%fist_env)) THEN
!          CALL fist_nonbond_env_get(f_env%force_env%fist_env%fist_nonbond_env, &
!                 rlist_lowsq=rlist_lowsq)
!        END IF
!      END IF
!    END IF
!    
!  END SUBROUTINE get_mimal_distances

END MODULE tmc_worker
