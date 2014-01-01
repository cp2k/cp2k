!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief module contains the master routine handling the tree creation,
!>        communication with workers and task distribution
!>        For each idle working group the master creates a new global tree 
!>        element, and if neccessay a related sub tree element, 
!>        OR find the next element to calculate the exact energy. 
!>        Goal is to keep at least the exact energy calculation working groups
!>        as busy as possible. 
!>        Master also checks for incomming results and update the tree and the
!>        acceptance ratios.
!> \par History
!>      11.2012 created [Mandes Schönherr]
!> \author Mandes
! *****************************************************************************

MODULE tmc_master
  USE cell_types,                      ONLY: init_cell
  USE cp_external_control,             ONLY: external_control
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush,&
                                             m_memory,&
                                             m_walltime
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE tmc_calculations,                ONLY: get_subtree_efficiency
  USE tmc_cancelation,                 ONLY: free_cancelation_list
  USE tmc_dot_tree,                    ONLY: create_dot_color,&
                                             create_global_tree_dot_color,&
                                             finalize_draw_tree,&
                                             init_draw_trees
  USE tmc_file_io,                     ONLY: print_restart_file,&
                                             write_element_in_file
  USE tmc_messages,                    ONLY: communicate_atom_types,&
                                             recv_msg,&
                                             send_msg,&
                                             stop_whole_group,&
                                             tmc_message
  USE tmc_move_handle,                 ONLY: check_moves,&
                                             print_move_types
  USE tmc_stati,                       ONLY: &
       TMC_CANCELING_MESSAGE, TMC_CANCELING_RECEIPT, TMC_STATUS_FAILED, &
       TMC_STATUS_WAIT_FOR_NEW_TASK, TMC_STATUS_WORKER_INIT, &
       TMC_STAT_ANALYSIS_REQUEST, TMC_STAT_ANALYSIS_RESULT, &
       TMC_STAT_APPROX_ENERGY_REQUEST, TMC_STAT_APPROX_ENERGY_RESULT, &
       TMC_STAT_ENERGY_REQUEST, TMC_STAT_ENERGY_RESULT, &
       TMC_STAT_INIT_ANALYSIS, TMC_STAT_MD_REQUEST, TMC_STAT_MD_RESULT, &
       TMC_STAT_NMC_REQUEST, TMC_STAT_NMC_RESULT, &
       TMC_STAT_SCF_STEP_ENER_RECEIVE, TMC_STAT_START_CONF_REQUEST, &
       TMC_STAT_START_CONF_RESULT
  USE tmc_tree_acceptance,             ONLY: check_acceptance_of_depending_subtree_nodes,&
                                             check_elements_for_acc_prob_update,&
                                             tree_update
  USE tmc_tree_build,                  ONLY: create_new_gt_tree_node,&
                                             deallocate_sub_tree_node,&
                                             finalize_init,&
                                             finalize_trees,&
                                             init_tree_mod,&
                                             remove_all_trees
  USE tmc_tree_search,                 ONLY: count_nodes_in_trees,&
                                             count_prepared_nodes_in_trees,&
                                             search_next_energy_calc
  USE tmc_tree_types,                  ONLY: &
       elem_array_type, elem_list_type, global_tree_type, status_accepted, &
       status_calc_approx_ener, status_calculate_MD, &
       status_calculate_NMC_steps, status_calculate_energy, &
       status_calculated, status_cancel_ener, status_cancel_nmc, &
       status_canceled_ener, status_canceled_nmc, status_created, &
       status_rejected, tree_type
  USE tmc_types,                       ONLY: tmc_env_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tmc_master'

  PUBLIC :: do_tmc_master

  INTEGER, PARAMETER :: DEBUG = 0

CONTAINS

! *****************************************************************************
!> \brief send cancel request to all workers processing elements in the list
!> \param cancel_list list with elements to cancel
!> \param work_list list with all elements processed by working groups
!> \param cancel_count counter of canceled elements
!> \param para_env communication environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE cancel_calculations(cancel_list, work_list, cancel_count,&
                                 para_env, tmc_env, error)
    TYPE(elem_list_type), POINTER            :: cancel_list
    TYPE(elem_array_type), DIMENSION(:), &
      POINTER                                :: work_list
    INTEGER                                  :: cancel_count
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(tmc_env_type), POINTER              :: tmc_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cancel_calculations', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat, wg
    LOGICAL                                  :: failure
    TYPE(elem_list_type), POINTER            :: tmp_element

    failure = .FALSE.
    IF(.NOT.ASSOCIATED(cancel_list)) RETURN
    NULLIFY(tmp_element)

    CPPrecondition(ASSOCIATED(tmc_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_env%params),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_env%m_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(work_list),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)


    stat = TMC_STATUS_FAILED
    wg   = -1
    cancel_elem_loop: DO
      ! find certain working group calculating this element
      working_elem_loop: DO i=1, SIZE(work_list)
        ! in special cases element could be distributed to several working groups, 
        ! but all, exept of one, should already be in canceling process
        IF((.NOT.work_list(i)%canceled).AND.&
           ASSOCIATED(work_list(i)%elem)) THEN
          IF(ASSOCIATED(cancel_list%elem, work_list(i)%elem)) THEN
            stat = TMC_CANCELING_MESSAGE
            wg = i
            EXIT working_elem_loop
          END IF
        END IF
      END DO working_elem_loop

      CPPostcondition(wg.GE.0,cp_failure_level,routineP,error,failure)
      CPPostcondition(stat.NE.TMC_STATUS_FAILED,cp_failure_level,routineP,error,failure)
      CPPostcondition(work_list(wg)%elem%stat.NE.status_calc_approx_ener,cp_failure_level,routineP,error,failure)

      
      IF(DEBUG.GE.1) &
        WRITE(tmc_env%m_env%io_unit,*)&
          "TMC|master: cancel group "//cp_to_string(wg)
      CALL tmc_message(msg_type=stat, send_recv=send_msg, dest=wg, &
               para_env=para_env, tmc_params=tmc_env%params, error=error)
      work_list(wg)%canceled = .TRUE.

      ! counting the amount of canceled elements
      cancel_count = cancel_count + 1

      ! delete element from canceling list
      IF(.NOT.ASSOCIATED(cancel_list%next)) THEN
         DEALLOCATE(cancel_list)
         cancel_list => NULL()
         EXIT cancel_elem_loop
      ELSE
         tmp_element => cancel_list%next
         DEALLOCATE(cancel_list)
         cancel_list => tmp_element
      END IF
    END DO cancel_elem_loop
  END SUBROUTINE cancel_calculations

! *****************************************************************************
!> \brief send analysis request to a worker 
!> \param ana_list list with elements to be analysed
!> \param work_list list with all elements processed by working groups
!> \param cancel_count counter of canceled elements
!> \param para_env communication environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE send_analysis_tasks(ana_list, ana_worker_info, para_env, tmc_env, error)
    TYPE(elem_list_type), POINTER            :: ana_list
    TYPE(elem_array_type), DIMENSION(:), &
      POINTER                                :: ana_worker_info
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(tmc_env_type), POINTER              :: tmc_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'send_analysis_tasks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat, wg
    LOGICAL                                  :: failure
    TYPE(elem_list_type), POINTER            :: list_tmp

    NULLIFY(list_tmp)

    CPPrecondition(ASSOCIATED(ana_worker_info),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)


    wg_loop:DO wg=1, SIZE(ana_worker_info)
      IF(.NOT.ASSOCIATED(ana_list)) EXIT wg_loop
      IF(.NOT.ana_worker_info(wg)%busy) THEN
        stat = TMC_STAT_ANALYSIS_REQUEST
        CALL tmc_message(msg_type=stat, send_recv=send_msg, dest=wg, &
               para_env=para_env, tmc_params=tmc_env%params, &
               list_elem=ana_list, error=error)
        IF(.NOT.ASSOCIATED(ana_list%next)) THEN
          DEALLOCATE(ana_list)
          ana_list => NULL()
        ELSE
          list_tmp => ana_list%next
          DEALLOCATE(ana_list)
          ana_list => list_tmp
        END IF
      END IF
    END DO wg_loop
  END SUBROUTINE send_analysis_tasks

! *****************************************************************************
!> \brief global master handling tree creation and communication/work 
!>        distribution with workers
!> \param tmc_env structure for storing all the tmc parameters
!> \param globenv global environment for external control
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE do_tmc_master(tmc_env, globenv, error)
    TYPE(tmc_env_type), POINTER              :: tmc_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'do_tmc_master', &
      routineP = moduleN//':'//routineN

    INTEGER :: cancel_count, handle, last_output, reactivation_cc_count, &
      reactivation_ener_count, restart_count, restarted_elem_nr, stat, &
      walltime_delay, walltime_offset, wg, worker_counter
    INTEGER, DIMENSION(6)                    :: nr_of_job
    INTEGER, DIMENSION(:), POINTER           :: tree_elem_counters, &
                                                tree_elem_heads
    LOGICAL                                  :: external_stop, failure, flag, &
                                                l_update_tree
    REAL(KIND=dp)                            :: run_time_start
    REAL(KIND=dp), DIMENSION(4)              :: worker_timings_aver
    REAL(KIND=dp), DIMENSION(:), POINTER     :: efficiency
    TYPE(elem_array_type), DIMENSION(:), &
      POINTER                                :: ana_worker_info, worker_info
    TYPE(global_tree_type), POINTER          :: gt_elem_tmp
    TYPE(tree_type), POINTER                 :: init_conf

    failure = .FALSE.
    external_stop = .FALSE.
    restarted_elem_nr = 0
    NULLIFY(init_conf, worker_info, ana_worker_info, gt_elem_tmp, tree_elem_counters)

    CPPrecondition(ASSOCIATED(tmc_env),cp_failure_level,routineP,error,failure)

    CPPrecondition(tmc_env%tmc_comp_set%group_nr==0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_env%tmc_comp_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_env%tmc_comp_set%para_env_m_w),cp_failure_level,routineP,error,failure)

    CPPrecondition(ASSOCIATED(tmc_env%m_env),cp_failure_level,routineP,error,failure)

    !-- run time measurment, to end just in time
    ! start the timing
    CALL timeset(routineN,handle)
    run_time_start = m_walltime()
    walltime_delay = 0
    walltime_offset = 20 ! default value the whole program needs to finalize

    ! initialize the different modules
    IF(tmc_env%params%DRAW_TREE) &
      CALL init_draw_trees(tmc_params=tmc_env%params, error=error) 

    !-- initialize variables
    ! nr_of_job: counting the differnt task send / received
    !  (1:NMC submitted, 2:energies submitted, 3:NMC finished 4:energy finished, 5:NMC canceled, 6:energy canceled)
    nr_of_job(:) = 0
    worker_counter = -1
    reactivation_ener_count = 0
    reactivation_cc_count = 0
    cancel_count = 0
    tmc_env%m_env%result_count = 0
    l_update_tree = .FALSE.
    restart_count = 1
    last_output = -1
    ! average timings
    !  (1:calculated NMC, 2:calculated ener, 3:canceled NMC, 4: canceled ener)
    worker_timings_aver(:) = 0.0_dp 
    ! remembers state of workers and their actual configurations
    ! the actual working group, communicating with
    ALLOCATE(worker_info(tmc_env%tmc_comp_set%para_env_m_w%num_pe-1))
    ALLOCATE(ana_worker_info(tmc_env%tmc_comp_set%para_env_m_ana%num_pe-1))

    ! get the start configuration form the first (exact energy) worker, 
    !   master should/could have no Force environment
    stat = TMC_STAT_START_CONF_REQUEST
    wg = 1
    CALL tmc_message(msg_type=stat, send_recv=send_msg, dest=wg, &
                     para_env=tmc_env%tmc_comp_set%para_env_m_w,&
                     tmc_params=tmc_env%params, &
                     wait_for_message=.TRUE., error=error)
    !-- wait for start configuration results and number of dimentions
    !-- get start configuration (init_conf element should not be allocated already)
    CALL tmc_message(msg_type=stat, send_recv=recv_msg, dest=wg, &
                     para_env=tmc_env%tmc_comp_set%para_env_m_w,&
                     tmc_params=tmc_env%params, &
                     elem=init_conf, success=flag, wait_for_message=.TRUE.,&
                     error=error)
    CALL cp_assert(stat.EQ.TMC_STAT_START_CONF_RESULT, cp_failure_level,&
                   cp_assertion_failed,routineP,&
                   "receiving start configuration failed, received stat "//&
                   cp_to_string(stat), error, failure=failure)
    IF(failure)THEN
      CALL stop_whole_group(para_env=tmc_env%tmc_comp_set%para_env_m_w, &
                            worker_info=worker_info, &
                            tmc_params=tmc_env%params, error=error)
      CALL stop_whole_group(para_env=tmc_env%tmc_comp_set%para_env_m_ana, &
                            worker_info=ana_worker_info, &
                            tmc_params=tmc_env%params, error=error)
      CALL timestop(handle)
      RETURN
    END IF
    ! get the atom names from first energy worker
    CALL communicate_atom_types(atoms=tmc_env%params%atoms, &
           source=1,&
           para_env=tmc_env%tmc_comp_set%para_env_m_first_w, &
           error=error)

    CALL init_cell(cell=tmc_env%params%cell)

    ! check the configuration consitency with selected moves
    CALL check_moves(tmc_params=tmc_env%params, &
           move_types=tmc_env%params%move_types, &
           mol_array=init_conf%mol, error=error)
    IF(ASSOCIATED(tmc_env%params%nmc_move_types)) &
      CALL check_moves(tmc_params=tmc_env%params, &
             move_types=tmc_env%params%nmc_move_types, &
             mol_array=init_conf%mol, error=error)

    ! set initial configuration 
    ! set initial random number generator seed (rng seed)
    ! initialize the tree structure espacially for parallel tmepering, 
    !   seting the subtrees
    CALL init_tree_mod(start_elem=init_conf, tmc_env=tmc_env, &
                       job_counts=nr_of_job, &
                       worker_timings=worker_timings_aver, error=error)

    ! init restart counter (espacially for restart case)
    IF(tmc_env%m_env%restart_out_step .NE. 0) THEN
      restart_count = INT(tmc_env%m_env%result_count(0)/&
                          REAL(tmc_env%m_env%restart_out_step, KIND=dp))+1
    END IF
    restarted_elem_nr = tmc_env%m_env%result_count(0)

!TODO check conf and cell of both input files (cell has to be equal, 
!           beacuse it is used as reference cell for scaling the cell)
    ! communicate the reference cell size
    DO wg=1, tmc_env%tmc_comp_set%para_env_m_w%num_pe-1
      stat = TMC_STATUS_WORKER_INIT
      CALL tmc_message(msg_type=stat, send_recv=send_msg, dest=wg, &
             para_env=tmc_env%tmc_comp_set%para_env_m_w,&
             tmc_params=tmc_env%params, &
             error=error)
    END DO    

    ! send the atom informations to all analysis workers
    IF(tmc_env%tmc_comp_set%para_env_m_ana%num_pe.GT.1)THEN
      DO wg=1, tmc_env%tmc_comp_set%para_env_m_ana%num_pe-1
        stat = TMC_STAT_INIT_ANALYSIS
        CALL tmc_message(msg_type=stat, send_recv=send_msg, dest=wg, &
               para_env=tmc_env%tmc_comp_set%para_env_m_ana,&
               result_count=tmc_env%m_env%result_count, &
               tmc_params=tmc_env%params, &
               elem=init_conf, &
               wait_for_message=.TRUE., error=error)
      END DO
      CALL communicate_atom_types(atoms=tmc_env%params%atoms, &
           source=0,&
           para_env=tmc_env%tmc_comp_set%para_env_m_ana, &
           error=error)
    END IF

    CALL deallocate_sub_tree_node(tmc_params=tmc_env%params, &
           tree_elem=init_conf, error=error)

    ! regtest output
    IF(tmc_env%params%print_test_output.OR.DEBUG.GT.0)&
      WRITE(tmc_env%m_env%io_unit,*)"TMC|first_global_tree_rnd_nr_X= ", &
                                    tmc_env%m_env%gt_head%rnd_nr

    ! calculate the approx energy of the first element (later the exact)
    IF(tmc_env%m_env%gt_head%conf(1)%elem%stat.EQ.status_calc_approx_ener)THEN
      wg = 1
      IF(tmc_env%tmc_comp_set%group_cc_nr.GT.0) &
        wg = tmc_env%tmc_comp_set%group_ener_nr +1
      stat = TMC_STAT_APPROX_ENERGY_REQUEST
      CALL tmc_message(msg_type=stat, send_recv=send_msg, dest=wg, &
                       para_env=tmc_env%tmc_comp_set%para_env_m_w,&
                       tmc_params=tmc_env%params, &
                       elem=tmc_env%m_env%gt_head%conf(1)%elem, &
                       error=error)
      worker_info(wg)%busy = .TRUE.
      worker_info(wg)%elem => tmc_env%m_env%gt_head%conf(1)%elem
      init_conf => tmc_env%m_env%gt_head%conf(1)%elem
    ELSE IF(tmc_env%m_env%gt_head%conf(1)%elem%stat.EQ.status_created)THEN
      init_conf => tmc_env%m_env%gt_head%conf(1)%elem
      ! calculation will be done automatically, 
      !   by searching the next conf for energy calculation
    END IF
    !-- START WORK --!
    !-- distributing work:
    !   1. receive incoming results
    !   2. check new results in tree
    !   3. if idle worker, create new tree element and send them to worker
    task_loop: DO
      ! =======================================================================
      !-- RECIEVING ALL incoming messages and handling them 
      ! results of tree node 1 is distributed to all other subtree nodes
      ! =======================================================================
      worker_request_loop: DO 
        wg = 1
        flag = .FALSE.
        CALL tmc_message(msg_type=stat,send_recv=recv_msg,dest=wg,&
                         para_env=tmc_env%tmc_comp_set%para_env_m_w,&
                         tmc_params=tmc_env%params, &
                         elem_array=worker_info(:), success=flag, error=error)

        IF(flag.EQV..FALSE.) EXIT worker_request_loop
        ! messages from worker group could be faster then the canceling request
        IF(worker_info(wg)%canceled .AND. (stat.NE.TMC_CANCELING_RECEIPT)) THEN
          IF(DEBUG.GE.1) &
            WRITE(tmc_env%m_env%io_unit,*)&
              "TMC|master: recv stat "//cp_to_string(stat)//&
              " of canceled worker group"
          CYCLE worker_request_loop
        END IF

        ! in case of parallel tempering canceled element could be reactivated,
        !   calculated faster and deleted
        CALL cp_assert(ASSOCIATED(worker_info(wg)%elem),&
                       cp_failure_level,cp_assertion_failed,routineP,&
                       "no tree elem exist when receiving stat "//&
                       cp_to_string(stat)//"of group"//cp_to_string(wg),&
                       error, failure=flag)

        IF(DEBUG.GE.1) &
          WRITE(tmc_env%m_env%io_unit,*)&
            "TMC|master: received stat "//cp_to_string(stat)// &
            " of sub tree "//cp_to_string(worker_info(wg)%elem%sub_tree_nr)//&
            " elem"//cp_to_string(worker_info(wg)%elem%nr)//&
            " with stat"//cp_to_string(worker_info(wg)%elem%stat)//&
            " of group"//cp_to_string(wg)//" group canceled ", worker_info(wg)%canceled
        SELECT CASE(stat)
        ! -- FAILED --------------------------
        CASE(TMC_STATUS_FAILED)
           EXIT task_loop 
        ! -- CANCEL_RECEIPT ------------------
        CASE(TMC_CANCELING_RECEIPT)
          ! worker should got cancel message before
          CPPrecondition(worker_info(wg)%canceled ,cp_failure_level,routineP,error,failure)
          worker_info(wg)%canceled = .FALSE.
          worker_info(wg)%busy = .FALSE.

          IF(ASSOCIATED(worker_info(wg)%elem)) THEN
            SELECT CASE(worker_info(wg)%elem%stat)
            CASE(status_cancel_ener)
               !-- timings
               worker_timings_aver(4) = (worker_timings_aver(4)*nr_of_job(6)+&
                 (m_walltime()-worker_info(wg)%start_time))/REAL(nr_of_job(6)+1,KIND=dp)
               nr_of_job(6) = nr_of_job(6) +1

               worker_info(wg)%elem%stat = status_canceled_ener
               worker_info(wg)%elem%potential = 8000.0_dp
               IF(tmc_env%params%DRAW_TREE) THEN
                 CALL create_dot_color(tree_element=worker_info(wg)%elem, &
                                       tmc_params=tmc_env%params, error=error)
               END IF
            CASE(status_cancel_nmc)
               !-- timings
               worker_timings_aver(3) = (worker_timings_aver(3)*nr_of_job(5)+&
                 (m_walltime()-worker_info(wg)%start_time))/REAL(nr_of_job(5)+1,KIND=dp)
               nr_of_job(5) = nr_of_job(5) +1

               worker_info(wg)%elem%stat = status_canceled_nmc
               worker_info(wg)%elem%potential = 8000.0_dp
               IF(tmc_env%params%DRAW_TREE) THEN 
                 CALL create_dot_color(tree_element=worker_info(wg)%elem, &
                                       tmc_params=tmc_env%params, error=error)
               END IF
            CASE DEFAULT 
               ! the subtree element is again in use (reactivated)
            END SELECT
            worker_info(wg)%elem => NULL()
          END IF
        ! -- START_CONF_RESULT ---------------
        CASE(TMC_STAT_START_CONF_RESULT)
           ! start configuration should already be handeled
           CPPrecondition(.FALSE. ,cp_failure_level,routineP,error,failure)
        ! -- ENERGY RESULT -----------------
        CASE(TMC_STAT_APPROX_ENERGY_RESULT)
          nr_of_job(3) = nr_of_job(3) +1
          worker_info(wg)%busy = .FALSE.
          worker_info(wg)%elem%stat = status_created
          IF(tmc_env%params%DRAW_TREE)THEN
            CALL create_dot_color(tree_element=worker_info(wg)%elem, &
                                  tmc_params=tmc_env%params, error=error)
          END IF
          worker_info(wg)%elem => NULL()
          ! nothing to do, the approximate potential 
          !   should be updated in the message interface 
        ! -- NMC / MD RESULT -----------------
        CASE(TMC_STAT_NMC_RESULT, TMC_STAT_MD_RESULT)
          IF(.NOT.worker_info(wg)%canceled) worker_info(wg)%busy = .FALSE.
          !-- timings for Nested Monte Carlo calculation
          worker_timings_aver(1) = (worker_timings_aver(1)*nr_of_job(3)+&
            (m_walltime()-worker_info(wg)%start_time))/REAL(nr_of_job(3)+1,KIND=dp)
          nr_of_job(3) = nr_of_job(3) +1

          worker_info(wg)%start_time = m_walltime() - worker_info(wg)%start_time
          CALL set_walltime_delay(worker_info(wg)%start_time, walltime_delay, error)
          worker_info(wg)%elem%stat = status_created
          IF(tmc_env%params%DRAW_TREE)THEN
            CALL create_dot_color(tree_element=worker_info(wg)%elem, &
                                  tmc_params=tmc_env%params, error=error)
          END IF
          !-- send energy request
          ! in case of one singe input file, energy is already calculated
          IF(tmc_env%params%NMC_inp_file.EQ."")THEN
             worker_info(wg)%elem%potential = worker_info(wg)%elem%e_pot_approx
             worker_info(wg)%elem%stat = status_calculated
             ! check acceptance of depending nodes
             IF(.NOT.(ASSOCIATED(worker_info(wg)%elem,init_conf))) THEN
               CALL check_acceptance_of_depending_subtree_nodes(tree_elem=worker_info(wg)%elem, &
                        tmc_env=tmc_env, error=error)
             END IF
             IF(tmc_env%params%DRAW_TREE)THEN
               CALL create_dot_color(tree_element=worker_info(wg)%elem, &
                                     tmc_params=tmc_env%params, error=error)
             END IF
             !-- CANCELING the calculations of the elements, which are definetively not needed anymore
             CALL cancel_calculations(cancel_list=tmc_env%m_env%cancelation_list, &
                                      work_list=worker_info, &
                                      para_env=tmc_env%tmc_comp_set%para_env_m_w, &
                                      tmc_env=tmc_env, &
                                      cancel_count=cancel_count, error=error)
             worker_info(wg)%elem => NULL()
          ELSE
             ! if all working groups are equal, the same group calculates the energy
             IF(tmc_env%tmc_comp_set%group_cc_nr.LE.0 &
                .AND.(.NOT.worker_info(wg)%canceled))THEN
               worker_info(wg)%elem%stat = status_calculate_energy
               stat = TMC_STAT_ENERGY_REQUEST
               ! immediately send energy request
               CALL tmc_message(msg_type=stat, send_recv=send_msg, dest=wg, &
                                para_env=tmc_env%tmc_comp_set%para_env_m_w, &
                                tmc_params=tmc_env%params, &
                                elem=worker_info(wg)%elem, error=error)
               worker_info(wg)%busy = .TRUE.
               nr_of_job(2) = nr_of_job(2) +1
               IF(tmc_env%params%DRAW_TREE)THEN
                 CALL create_dot_color(tree_element=worker_info(wg)%elem, &
                                       tmc_params=tmc_env%params, error=error)
               END IF
               !-- set start time for energy calculation
               worker_info(wg)%start_time = m_walltime() 
             ELSE 
               worker_info(wg)%elem => NULL()
             END IF
          END IF
        ! -- ENERGY RESULT --------------------
        CASE(TMC_STAT_ENERGY_RESULT)
          !-- timings
          worker_timings_aver(2) = (worker_timings_aver(2)*nr_of_job(4)+&
            (m_walltime()-worker_info(wg)%start_time))/REAL(nr_of_job(4)+1,KIND=dp)
          nr_of_job(4) = nr_of_job(4) +1

          worker_info(wg)%start_time = m_walltime()-worker_info(wg)%start_time
          CALL set_walltime_delay(worker_info(wg)%start_time, walltime_delay, error)

          IF(.NOT.worker_info(wg)%canceled)&
            worker_info(wg)%busy = .FALSE.
          ! the first node in tree is always accepted.!.
          IF(ASSOCIATED(worker_info(wg)%elem, init_conf)) THEN
            !-- distribute energy of first element to all subtrees
            CALL finalize_init(gt_tree_ptr=tmc_env%m_env%gt_head, &
                               tmc_env=tmc_env, error=error)
            IF(tmc_env%params%DRAW_TREE) THEN
              CALL create_global_tree_dot_color(gt_tree_element=tmc_env%m_env%gt_act, &
                                                tmc_params=tmc_env%params, error=error)
              CALL create_dot_color(tree_element=worker_info(wg)%elem, &
                                    tmc_params=tmc_env%params, error=error)
            END IF
            init_conf => NULL()
          ELSE
            worker_info(wg)%elem%stat = status_calculated
            IF(tmc_env%params%DRAW_TREE) &
              CALL create_dot_color(worker_info(wg)%elem, &
                                    tmc_params=tmc_env%params, error=error)
            ! check acceptance of depending nodes
            ! first (initial) configuration do not have to be checked
            CALL check_acceptance_of_depending_subtree_nodes(tree_elem=worker_info(wg)%elem, &
                     tmc_env=tmc_env,error=error)
          END IF
          !-- write out all configurations (not only Markov Chain) e.g. for fitting
          IF(tmc_env%params%all_conf_file_name.NE."")THEN
            CALL write_element_in_file(elem=worker_info(wg)%elem, &
                                       file_name=tmc_env%params%all_conf_file_name, &
                                       tmc_params=tmc_env%params, &
                                       conf_nr=nr_of_job(4), error=error)
          END IF

          !-- CANCELING the calculations of the elements, 
          !      which are definetively not needed anymore
          CALL cancel_calculations(cancel_list=tmc_env%m_env%cancelation_list,&
                                   work_list=worker_info, &
                                   para_env=tmc_env%tmc_comp_set%para_env_m_w, &
                                   tmc_env=tmc_env, &
                                   cancel_count=cancel_count, error=error)
          IF(DEBUG.GE.9) &
            WRITE(tmc_env%m_env%io_unit,*) &
              "TMC|master: handled energy result of sub tree ", &
              worker_info(wg)%elem%sub_tree_nr , " elem ", worker_info(wg)%elem%nr,&
              " with stat",  worker_info(wg)%elem%stat
          worker_info(wg)%elem => NULL()

        !-- SCF ENERGY -----------------------
        CASE(TMC_STAT_SCF_STEP_ENER_RECEIVE)
          IF(.NOT.(ASSOCIATED(worker_info(wg)%elem,init_conf)).AND.&
            worker_info(wg)%elem%stat.NE.status_cancel_ener.AND.&
            worker_info(wg)%elem%stat.NE.status_cancel_nmc) THEN
            ! update the acceptance probability and the canceling list 
            CALL check_elements_for_acc_prob_update(tree_elem=worker_info(wg)%elem, &
                                                    tmc_env=tmc_env, error=error)
          END IF
          ! cancel inlikely elements
          CALL cancel_calculations(cancel_list=tmc_env%m_env%cancelation_list, &
                                   work_list=worker_info, &
                                   para_env=tmc_env%tmc_comp_set%para_env_m_w, &
                                   tmc_env=tmc_env, &
                                   cancel_count=cancel_count,error=error)
        CASE(TMC_STAT_ANALYSIS_RESULT)
          ana_worker_info(wg)%busy = .FALSE.
          ana_worker_info(wg)%elem => NULL()
        CASE DEFAULT
          CALL cp_assert(.FALSE.,&
                   cp_failure_level,cp_assertion_failed,routineP,&
                   "received message with unknown info/stat type", &
                   error, failure=flag)
        END SELECT
      END DO worker_request_loop
      !-- do tree update (check new results)
      CALL tree_update(tmc_env=tmc_env, result_acc=flag, &
                       something_updated=l_update_tree, error=error)
      IF(DEBUG.GE.2 .AND. l_update_tree) &
        WRITE(tmc_env%m_env%io_unit,*) &
          "TMC|master: tree updated "//cp_to_string(l_update_tree)// &
          " of with gt elem "//cp_to_string(tmc_env%m_env%gt_act%nr)//&
          " with stat"//cp_to_string(tmc_env%m_env%gt_act%stat)

      CALL send_analysis_tasks(ana_list=tmc_env%m_env%analysis_list, &
                               ana_worker_info=ana_worker_info, &
                               para_env=tmc_env%tmc_comp_set%para_env_m_ana, &
                               tmc_env=tmc_env, error=error)

      ! =======================================================================
      !-- ALL CALCULATIONS DONE (check) ---
      ! =======================================================================
      ! if enough configurations are sampled or walltime is exeeded, 
      !   finish building trees
!TODO set correct logger para_env to use this 
      CALL external_control(should_stop=external_stop,flag="TMC",globenv=globenv,&
                            error=error)
      IF((ANY(tmc_env%m_env%result_count(1:).GE.tmc_env%m_env%num_MC_elem)&
         .AND.flag).OR.&
         (m_walltime()-run_time_start .GT. &
         tmc_env%m_env%walltime-walltime_delay-walltime_offset).OR. &
         external_stop) THEN
        WRITE(tmc_env%m_env%io_unit,FMT="(/,T2,A)") REPEAT("=",79)
        ! calculations NOT finished, walltime exceeded
        IF(.NOT.ANY(tmc_env%m_env%result_count(1:) &
                    .GE.tmc_env%m_env%num_MC_elem)) THEN
          WRITE(tmc_env%m_env%io_unit,*)"Walltime exceeded.", &
                m_walltime()-run_time_start," of ", tmc_env%m_env%walltime-walltime_delay-walltime_offset, &
                "(incl. delay", walltime_delay, "and offset", walltime_offset,") left"
        ELSE
        ! calculations finished
          IF(tmc_env%params%print_test_output)&
            WRITE(tmc_env%m_env%io_unit,*)"Total energy: ", &
              tmc_env%m_env%result_list(1)%elem%potential
        END IF
        IF(tmc_env%m_env%restart_out_step.NE.0)&
          CALL print_restart_file(tmc_env=tmc_env,job_counts=nr_of_job, &
                                  timings=worker_timings_aver, error=error)
        EXIT task_loop
      END IF

      ! =======================================================================
      ! update the rest of the tree (canceling and deleting elements)
      ! =======================================================================
      IF(l_update_tree) THEN
        IF(DEBUG.GE.2) &
          WRITE(tmc_env%m_env%io_unit,*)&
            "TMC|master: start remove elem and cancel calculation"
        !-- CLEANING tree nodes beside the path through the tree from 
        !      end_of_clean_tree to tree_ptr
        ! --> getting back the end of clean tree
        CALL remove_all_trees(working_elem_list=worker_info, tmc_env=tmc_env, &
                               error=error) 
        !-- CANCELING the calculations of the elements, 
        !     which are definetively not needed anymore
        !   elements are added to canceling list if no global tree reference 
        !     exist anymore
        CALL cancel_calculations(cancel_list=tmc_env%m_env%cancelation_list, &
                                 work_list=worker_info, &
                                 cancel_count=cancel_count, &
                                 para_env=tmc_env%tmc_comp_set%para_env_m_w,&
                                 tmc_env=tmc_env, error=error)
      END IF

      ! =====================================================================
      !-- NEW TASK (if worker not busy sumit next task)
      ! =====================================================================
      worker_counter = worker_counter + 1
      wg = MODULO(worker_counter, tmc_env%tmc_comp_set%para_env_m_w%num_pe-1)+1

      IF(DEBUG.GE.16.AND.ALL(worker_info(:)%busy))&
        WRITE(tmc_env%m_env%io_unit,*) "all workers are busy"

      IF(.NOT.worker_info(wg)%busy)THEN
        IF(DEBUG.GE.13) &
          WRITE(tmc_env%m_env%io_unit,*)&
            "TMC|master: search new task for worker ", wg
        ! no group separation
        IF(tmc_env%tmc_comp_set%group_cc_nr.LE.0)THEN
          ! search next element to calculate the energy
          CALL search_next_energy_calc(gt_head=tmc_env%m_env%gt_act, &
                                       new_gt_elem=gt_elem_tmp, stat=stat,&
                                       react_count=reactivation_ener_count, &
                                       error=error)
          IF(stat.EQ.TMC_STATUS_WAIT_FOR_NEW_TASK) THEN
            CALL create_new_gt_tree_node(tmc_env=tmc_env, stat=stat, &
                   new_elem=gt_elem_tmp, &
                   reactivation_cc_count=reactivation_cc_count, &
                   error=error)
          END IF
        ELSEIF(wg.GT.tmc_env%tmc_comp_set%group_ener_nr) THEN
        ! specialized groups (groups for exact energy and groups for configurational change)
          ! creating new element (configurational change group)
          !-- crate new node, configurational change is handled in tmc_tree module
          CALL create_new_gt_tree_node(tmc_env=tmc_env, stat=stat, &
                   new_elem=gt_elem_tmp, &
                   reactivation_cc_count=reactivation_cc_count, &
                   error=error)
          ! element could be already created, hence CC worker has nothing to do for this element
          ! in next round he will get a task
          IF(stat.EQ.status_created .OR. stat.EQ.status_calculate_energy) &
            stat = TMC_STATUS_WAIT_FOR_NEW_TASK
        ELSE
          ! search next element to calculate the energy
          CALL search_next_energy_calc(gt_head=tmc_env%m_env%gt_act, &
                                       new_gt_elem=gt_elem_tmp, stat=stat, &
                                       react_count=reactivation_ener_count, &
                                       error=error)
        END IF

        IF(DEBUG.GE.10)&
          WRITE(tmc_env%m_env%io_unit,*) &
            "TMC|master: send task with elem stat "//cp_to_string(stat)// &
            " to group "//cp_to_string(wg)
        ! MESSAGE settings: status informations and task for communication
        SELECT CASE(stat)
        CASE(TMC_STATUS_WAIT_FOR_NEW_TASK)
          CYCLE task_loop
        CASE(TMC_STATUS_FAILED)
          !STOP "in creating new task, status failed should be handled before"
          CYCLE task_loop
        CASE(status_calculated, status_accepted, status_rejected)
          CYCLE task_loop
        CASE(status_calc_approx_ener)
          ! e.g. after volume move, we need the approximate potential for 2 potential check of following NMC nodes
          stat = TMC_STAT_APPROX_ENERGY_REQUEST
          CALL tmc_message(msg_type=stat, send_recv=send_msg, dest=wg, &
                           para_env=tmc_env%tmc_comp_set%para_env_m_w,&
                           tmc_params=tmc_env%params, &
                           elem=gt_elem_tmp%conf(gt_elem_tmp%mv_conf)%elem, &
                           error=error)
          nr_of_job(1) = nr_of_job(1) +1
        CASE(status_created, status_calculate_energy)
          ! in case of parallel tempering the node can be already be calculating (related to another global tree node
          !-- send task to calculate system property
          gt_elem_tmp%conf(gt_elem_tmp%mv_conf)%elem%stat = status_calculate_energy
          IF(tmc_env%params%DRAW_TREE) &
            CALL create_dot_color(tree_element=gt_elem_tmp%conf(gt_elem_tmp%mv_conf)%elem, &
                                  tmc_params=tmc_env%params, error=error)
          stat = TMC_STAT_ENERGY_REQUEST
          CALL tmc_message(msg_type=stat, send_recv=send_msg, dest=wg, &
                           para_env=tmc_env%tmc_comp_set%para_env_m_w,&
                           tmc_params=tmc_env%params, &
                           elem=gt_elem_tmp%conf(gt_elem_tmp%mv_conf)%elem, &
                           error=error)
          nr_of_job(2) = nr_of_job(2) +1
        CASE(status_calculate_MD)
          stat = TMC_STAT_MD_REQUEST
          CALL tmc_message(msg_type=stat, send_recv=send_msg, dest=wg, &
                           para_env=tmc_env%tmc_comp_set%para_env_m_w,&
                           tmc_params=tmc_env%params, &
                           elem=gt_elem_tmp%conf(gt_elem_tmp%mv_conf)%elem, &
!                           temperature=tmc_env%params%Temp(gt_elem_tmp%mv_conf), &
                           error=error)
          nr_of_job(1) = nr_of_job(1) +1
        CASE(status_calculate_NMC_steps)
          !-- send information of element, which should be calculated
          stat = TMC_STAT_NMC_REQUEST
          CALL tmc_message(msg_type=stat, send_recv=send_msg, dest=wg, &
                           para_env=tmc_env%tmc_comp_set%para_env_m_w,&
                           tmc_params=tmc_env%params, &
                           elem=gt_elem_tmp%conf(gt_elem_tmp%mv_conf)%elem, &
                           error=error)
          nr_of_job(1) = nr_of_job(1) +1
        CASE(status_cancel_nmc, status_cancel_ener)
          ! skip that task until receipt is received
          ! no status update
        CASE DEFAULT
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                         routineP,"new task of tree element"//&
                         cp_to_string(gt_elem_tmp%nr)//&
                         "has unknown status"//cp_to_string(stat), &
                         error, failure=flag)
        END SELECT
        worker_info(wg)%elem => gt_elem_tmp%conf(gt_elem_tmp%mv_conf)%elem 
        worker_info(wg)%busy = .TRUE.
        ! set timer for maximum calculation time recognition 
        worker_info(wg)%start_time = m_walltime()

        !===================== write out info after x requested tasks==========
        IF(nr_of_job(4).GT.last_output .AND. &
           (MODULO(nr_of_job(4),tmc_env%m_env%info_out_step_size).EQ.0).AND.&
           (stat.NE.TMC_STATUS_FAILED)) THEN
          last_output = nr_of_job(4)
          WRITE(tmc_env%m_env%io_unit,FMT="(/,T2,A)") REPEAT("-",79)
          WRITE(tmc_env%m_env%io_unit,*) &
                  "Tasks submitted:  E ",nr_of_job(2), ", cc",nr_of_job(1)
          WRITE(tmc_env%m_env%io_unit,*) &
                  "Results received: E ",nr_of_job(4), ", cc",nr_of_job(3)
          WRITE(tmc_env%m_env%io_unit,*) &
                  "Configurations used:",tmc_env%m_env%result_count(0), &
                  ", sub trees", tmc_env%m_env%result_count(1:)

          CALL print_move_types(init=.FALSE., file_io=tmc_env%m_env%io_unit, &
                                tmc_params=tmc_env%params, error=error)
          ALLOCATE(tree_elem_counters(0:SIZE(tmc_env%params%Temp)))
          ALLOCATE(tree_elem_heads(0:SIZE(tmc_env%params%Temp)))
          CALL count_nodes_in_trees(global_tree_ptr=tmc_env%m_env%gt_act, &
                 end_of_clean_trees=tmc_env%m_env%st_clean_ends,&
                 counters=tree_elem_counters, head_elements_nr=tree_elem_heads, &
                 error=error) 
          WRITE(tmc_env%m_env%io_unit,*) "nodes in trees", tree_elem_counters(:)
          WRITE(tmc_env%m_env%io_unit,*) "tree heads    ", tree_elem_heads(:) 
          IF(tmc_env%params%NMC_inp_file.NE."") THEN
            CALL count_prepared_nodes_in_trees(global_tree_ptr=tmc_env%m_env%gt_act, &
                   counters=tree_elem_counters, error=error)
            WRITE(tmc_env%m_env%io_unit,FMT=*)&
                    "ener prepared ", tree_elem_counters
          END IF
          IF(tmc_env%params%SPECULATIVE_CANCELING) &
            WRITE(tmc_env%m_env%io_unit,*)&
                  "canceled cc|E:     ", nr_of_job(5:6), &
                  ", reactivated: cc ",&
                  reactivation_cc_count, &
                  ", reactivated: E ",  &
                  reactivation_ener_count
            WRITE(tmc_env%m_env%io_unit,FMT='(A,2F10.2)')&
                  " Average time for cc/ener calc          ", &
                  worker_timings_aver(1), worker_timings_aver(2)
          IF(tmc_env%params%SPECULATIVE_CANCELING) &
            WRITE(tmc_env%m_env%io_unit,FMT='(A,2F10.2)')&
                  " Average time until cancel cc/ener calc ", &
                  worker_timings_aver(3), worker_timings_aver(4)
          IF(tmc_env%params%esimate_acc_prob) &
            WRITE(tmc_env%m_env%io_unit,*)&
                  "Estimate correct (acc/Nacc) | wrong (acc/nacc)", &
                  tmc_env%m_env%estim_corr_wrong(1),&
                  tmc_env%m_env%estim_corr_wrong(3), " | ",&
                  tmc_env%m_env%estim_corr_wrong(2),&
                  tmc_env%m_env%estim_corr_wrong(4)
          WRITE(tmc_env%m_env%io_unit,*)&
                  "Time: ",INT(m_walltime()-run_time_start), "of",&
                  INT(tmc_env%m_env%walltime-walltime_delay-walltime_offset),&
                  "sec needed. "
          WRITE(tmc_env%m_env%io_unit,*)&
                  "Memory used: ", INT(m_memory()/(1024*1024),KIND=KIND(0)), "Mbytes"
          CALL m_flush(tmc_env%m_env%io_unit)
          DEALLOCATE(tree_elem_heads)
          DEALLOCATE(tree_elem_counters)
        END IF
        !===================== write out restart file after x results============
        IF(tmc_env%m_env%restart_out_step.GT.0 .AND. &
           tmc_env%m_env%result_count(0).GT.&
           restart_count*tmc_env%m_env%restart_out_step) THEN
          CALL print_restart_file(tmc_env=tmc_env, job_counts=nr_of_job, &
                                  timings=worker_timings_aver, error=error)
          restart_count = restart_count +1
        END IF

      END IF !worker busy?
    END DO task_loop

    ! -- END OF WORK (enough configurations are calculated or walltime exceeded
    WRITE(tmc_env%m_env%io_unit,FMT="(/,T2,A)") REPEAT("=",79)
    WRITE (UNIT=tmc_env%m_env%io_unit,FMT="(T2,A,T35,A,T80,A)") "=",&
      "finalizing TMC","="
    WRITE(tmc_env%m_env%io_unit,*) "acceptance rates:"
    CALL print_move_types(init=.FALSE., file_io=tmc_env%m_env%io_unit, &
                          tmc_params=tmc_env%params, error=error)
    WRITE(tmc_env%m_env%io_unit,FMT="(/,T2,A)") REPEAT("-",79)
    ! program efficency result outputs
    ALLOCATE(efficiency(0:tmc_env%params%nr_temp))
    efficiency = 0.0_dp
    CALL get_subtree_efficiency(tmc_env=tmc_env, eff=efficiency, error=error)
    WRITE(tmc_env%m_env%io_unit,*) "Efficiencies:" 
    WRITE(tmc_env%m_env%io_unit,FMT="(A,F5.2,A,1000F5.2)")&
      " (MC elements/calculated configuration) global:", &
      efficiency(0), " sub tree(s): ",efficiency(1:)
    DEALLOCATE(efficiency)
    IF(tmc_env%tmc_comp_set%group_cc_nr.GT.0) &
    WRITE(tmc_env%m_env%io_unit,FMT="(A,1000F5.2)")&
       " (MC elements/created configuration)          :", &
       tmc_env%m_env%result_count(:)/REAL(nr_of_job(3),KIND=dp)
    WRITE(tmc_env%m_env%io_unit,FMT="(A,1000F5.2)") &
      " (MC elements/energy calculated configuration):", &
      tmc_env%m_env%result_count(:)/REAL(nr_of_job(4),KIND=dp)
    IF(tmc_env%params%NMC_inp_file.NE."") THEN
      WRITE(tmc_env%m_env%io_unit,*) &
        "Amount of canceled elements (E/cc):",&
        tmc_env%m_env%count_cancel_ener, tmc_env%m_env%count_cancel_NMC
      WRITE(tmc_env%m_env%io_unit,*) &
        " reactivated E  ", reactivation_ener_count
      WRITE(tmc_env%m_env%io_unit,*) &
        " reactivated cc ", reactivation_cc_count
    END IF
    WRITE(tmc_env%m_env%io_unit,FMT="(A,F10.2)") &
      " computing time of one Markov chain element ", &
      (m_walltime()-run_time_start)/REAL(tmc_env%m_env%result_count(0)-&
                                         restarted_elem_nr,KIND=dp)
    WRITE(tmc_env%m_env%io_unit,FMT="(A,F10.2)") " TMC run time[s]: ", m_walltime()-run_time_start
    WRITE(tmc_env%m_env%io_unit,FMT="(/,T2,A)") REPEAT("=",79)

    !-- FINALIZE 
    WRITE(tmc_env%m_env%io_unit,*) "stopping workers"
    CALL stop_whole_group(para_env=tmc_env%tmc_comp_set%para_env_m_w, &
                          worker_info=worker_info, &
                          tmc_params=tmc_env%params, error=error)
    DEALLOCATE(worker_info)
    CALL stop_whole_group(para_env=tmc_env%tmc_comp_set%para_env_m_ana, &
                          worker_info=ana_worker_info, &
                          tmc_params=tmc_env%params, error=error)
    DEALLOCATE(ana_worker_info)

    !-- deallocating everything in tree module
    CALL finalize_trees(tmc_env=tmc_env, error=error)

    CALL free_cancelation_list(tmc_env%m_env%cancelation_list, error)

    ! -- write final configuration
    IF(tmc_env%params%DRAW_TREE) &
      CALL finalize_draw_tree(tmc_params=tmc_env%params, error=error)

    WRITE(tmc_env%m_env%io_unit,*) "TMC master: all work done."

    ! end the timing
    CALL timestop(handle) 

  END SUBROUTINE do_tmc_master

! *****************************************************************************
!> \brief routine sets the walltime delay, to the maximum calculation time
!>        hence the program can stop with a proper finailze
!> \param time actual calculation time
!> \param walltime_delay the actual biggest calculation time
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE set_walltime_delay(time, walltime_delay, error)
    REAL(KIND=dp)                            :: time
    INTEGER                                  :: walltime_delay
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_walltime_delay', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(time.GE.0.0_dp,cp_failure_level,routineP,error,failure)

    IF(time.GT.walltime_delay)THEN
      walltime_delay = INT(time)+1
    END IF
  END SUBROUTINE set_walltime_delay

END MODULE tmc_master
