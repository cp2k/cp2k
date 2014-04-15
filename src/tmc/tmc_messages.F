!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief set up the different message for different tasks
!>      A TMC message consists of 3 parts (messages)
!>      1: first a message with task type (STATUS) and SIZES of submessages
!>      2: (if existing) a message with INTEGER values
!>      3: (if existing) a message with REAL values
!>      submessages 2 and 3 include relevant data, e.g. positions, box sizes...
!> \par History
!>      11.2012 created [Mandes Schönherr]
!> \author Mandes
! *****************************************************************************
MODULE tmc_messages
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE message_passing,                 ONLY: mp_any_source,&
                                             mp_any_tag,&
                                             mp_bcast,&
                                             mp_probe,&
                                             mp_recv,&
                                             mp_send
  USE tmc_move_handle,                 ONLY: add_mv_prob
  USE tmc_stati,                       ONLY: &
       TMC_CANCELING_MESSAGE, TMC_CANCELING_RECEIPT, TMC_STATUS_CALCULATING, &
       TMC_STATUS_FAILED, TMC_STATUS_STOP_RECEIPT, &
       TMC_STATUS_WAIT_FOR_NEW_TASK, TMC_STATUS_WORKER_INIT, &
       TMC_STAT_ANALYSIS_REQUEST, TMC_STAT_ANALYSIS_RESULT, &
       TMC_STAT_APPROX_ENERGY_REQUEST, TMC_STAT_APPROX_ENERGY_RESULT, &
       TMC_STAT_ENERGY_REQUEST, TMC_STAT_ENERGY_RESULT, &
       TMC_STAT_INIT_ANALYSIS, TMC_STAT_MD_BROADCAST, TMC_STAT_MD_REQUEST, &
       TMC_STAT_MD_RESULT, TMC_STAT_NMC_BROADCAST, TMC_STAT_NMC_REQUEST, &
       TMC_STAT_NMC_RESULT, TMC_STAT_SCF_STEP_ENER_RECEIVE, &
       TMC_STAT_START_CONF_REQUEST, TMC_STAT_START_CONF_RESULT, &
       task_type_gaussian_adaptation
  USE tmc_tree_build,                  ONLY: allocate_new_sub_tree_node
  USE tmc_tree_types,                  ONLY: elem_array_type,&
                                             elem_list_type,&
                                             tree_type
  USE tmc_types,                       ONLY: allocate_tmc_atom_type,&
                                             tmc_atom_type,&
                                             tmc_param_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PARAMETER, PUBLIC                 :: send_msg = .TRUE.
  LOGICAL, PARAMETER, PUBLIC                 :: recv_msg = .FALSE.

  INTEGER, PARAMETER                         :: message_end_flag = 25

  INTEGER, PARAMETER                         :: DEBUG = 0
  INTEGER, PARAMETER                         :: tmc_out_file_nr =6

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tmc_messages'

  PUBLIC :: check_if_group_master
  PUBLIC :: tmc_message
  PUBLIC :: communicate_atom_types
  PUBLIC :: stop_whole_group

  INTEGER, PARAMETER, PUBLIC :: MASTER_COMM_ID = 0 ! id for master and group master
  INTEGER, PARAMETER, PUBLIC :: bcast_group    = -1! destination flag for broadcasting to other group participants
  INTEGER, PARAMETER :: TMC_SEND_INFO_SIZE     = 4 ! usually: 1. status, array sizes: 2. int, 3. real, 4. char

  TYPE message_send
    INTEGER, DIMENSION(TMC_SEND_INFO_SIZE)   :: info
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: task_real
    INTEGER, DIMENSION(:), ALLOCATABLE       :: task_int
    CHARACTER, DIMENSION(:), ALLOCATABLE     :: task_char
    !should be deleted somewhen
    INTEGER, DIMENSION(:), ALLOCATABLE :: elem_stat
  END TYPE message_send

CONTAINS

! *****************************************************************************
!> \brief checks if the core is the group master
!> \param para_env defines the mpi communicator
!> \param master return value, logical
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 01.2013
! *****************************************************************************
  FUNCTION check_if_group_master(para_env, error) RESULT(master)
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: master

    CHARACTER(LEN=*), PARAMETER :: routineN = 'check_if_group_master', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error, failure)

    master = .FALSE.
    IF(para_env%mepos.EQ.MASTER_COMM_ID) &
      master = .TRUE.
  END FUNCTION check_if_group_master

! *****************************************************************************
!> \brief tmc message handling, packing messages with integer and real data 
!>        type. Send first info message with task type and message sizes and
!>        then the int and real messages. The same for receiving
!> \param msg_type defines the message types, see message tags definition
!> \param send_recv: 1= send, 0= receive
!> \param dest: defines the target or source of message 
!>              (-1=braodcast, 0= master, 1... working group)
!> \param para_env defines the mpi communicator
!> \param tmc_params stuct with parameters (global settings)
!> \param elem: a subtree element from which info are readed or written in
!> \param T_ind: (optinal) given because not in tree structure
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE tmc_message(msg_type, send_recv, dest, para_env, tmc_params, &
                         elem, elem_array, list_elem, result_count, &
                         wait_for_message, success, error)
    INTEGER                                  :: msg_type
    LOGICAL                                  :: send_recv
    INTEGER                                  :: dest
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(tree_type), OPTIONAL, POINTER       :: elem
    TYPE(elem_array_type), DIMENSION(:), &
      OPTIONAL                               :: elem_array
    TYPE(elem_list_type), OPTIONAL, POINTER  :: list_elem
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: result_count
    LOGICAL, OPTIONAL                        :: wait_for_message, success
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tmc_message', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, message_tag, tmp_tag
    LOGICAL                                  :: act_send_recv, failure, flag
    TYPE(message_send), POINTER              :: m_send

    CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

    ALLOCATE(m_send)

    ! init
    ! define send_recv flag for broadcast
    IF(dest.EQ.bcast_group) THEN
      ! master should always send
      IF(para_env%mepos.EQ.MASTER_COMM_ID) THEN
        act_send_recv = send_msg
      ELSE
      ! worker should always receive
        act_send_recv = recv_msg
      END IF
    ELSE
      act_send_recv = send_recv
    END IF
    message_tag = 0

    ! =============================
    ! sending message
    ! =============================
    ! creating message to send
    IF(act_send_recv.EQV.send_msg) THEN
      IF((DEBUG.GE.7).AND.(dest.NE.bcast_group).AND.&
         (dest.NE.MASTER_COMM_ID)) THEN
         IF(PRESENT(elem)) THEN
           WRITE(tmc_out_file_nr,*)"send element info to ", dest, " of type ", msg_type,"of subtree", elem%sub_tree_nr, &
                                   "elem", elem%nr
         ELSE
           WRITE(tmc_out_file_nr,*)"send element info to ", dest, " of type ", msg_type
         END IF
      END IF
      SELECT CASE (msg_type)
      CASE(TMC_STAT_START_CONF_REQUEST, TMC_STATUS_FAILED, TMC_CANCELING_MESSAGE, &
           TMC_CANCELING_RECEIPT, TMC_STATUS_STOP_RECEIPT, &
           TMC_STATUS_WAIT_FOR_NEW_TASK, TMC_STATUS_CALCULATING, &
           TMC_STAT_ANALYSIS_RESULT)
        CALL create_status_message(msg_type, m_send, error)
      CASE(TMC_STATUS_WORKER_INIT)
        CALL create_worker_init_message(msg_type, tmc_params, m_send, error)
      CASE(TMC_STAT_START_CONF_RESULT, TMC_STAT_INIT_ANALYSIS)
        CALL create_start_conf_message(msg_type, elem, result_count, tmc_params, m_send, error)
      CASE(TMC_STAT_ENERGY_REQUEST, TMC_STAT_APPROX_ENERGY_REQUEST)
        CALL create_energy_request_message(msg_type, elem, m_send, tmc_params, error)
      CASE(TMC_STAT_APPROX_ENERGY_RESULT)
        CALL create_approx_energy_result_message(msg_type, elem, m_send, tmc_params, error)
      CASE(TMC_STAT_ENERGY_RESULT)
        CALL create_energy_result_message(msg_type, elem, m_send, tmc_params, error)
      CASE(TMC_STAT_NMC_REQUEST, TMC_STAT_NMC_BROADCAST, &
           TMC_STAT_MD_REQUEST, TMC_STAT_MD_BROADCAST)
        CALL create_NMC_request_massage(msg_type, elem, m_send, tmc_params, error)
      CASE(TMC_STAT_MD_RESULT, TMC_STAT_NMC_RESULT)
        CALL create_NMC_result_massage(msg_type, elem, m_send, tmc_params, error)
      CASE(TMC_STAT_ANALYSIS_REQUEST)
        CPPrecondition(PRESENT(list_elem),cp_failure_level,routineP,error,failure)
        CALL create_analysis_request_message(list_elem, m_send, tmc_params, error)
      CASE DEFAULT
        CALL cp_assert(.FALSE., cp_failure_level,cp_assertion_failed,routineP,&
                       "try to send unknown message type "//cp_to_string(msg_type),&
                       error, failure=flag)
      END SELECT
      !set message info
      message_tag = msg_type
      m_send%info(:) = 0
      m_send%info(1) = msg_type
      IF(ALLOCATED(m_send%task_int))  m_send%info(2) = SIZE(m_send%task_int)
      IF(ALLOCATED(m_send%task_real)) m_send%info(3) = SIZE(m_send%task_real)
      IF(ALLOCATED(m_send%task_char)) m_send%info(4) = SIZE(m_send%task_char)
    END IF

    ! sending message
    IF((act_send_recv.EQV.send_msg).AND.(dest.NE.bcast_group)) THEN
      CALL mp_send(m_send%info, dest,message_tag, para_env%group)
      IF(m_send%info(2).GT.0) THEN
        CALL mp_send(m_send%task_int, dest, message_tag, para_env%group)      
      END IF
      IF(m_send%info(3).GT.0) THEN
        CALL mp_send(m_send%task_real, dest, message_tag, para_env%group)
      END IF
      IF(m_send%info(4).GT.0) THEN
        CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
        !TODO send characters CALL mp_send(m_send%task_char, dest, message_tag, para_env%group)
      END IF
      IF(DEBUG.GE.1) &
        WRITE(tmc_out_file_nr,*)"TMC|message: ID: ",para_env%mepos, &
          " send element info to   ", dest, " of stat ", m_send%info(1), &
          " with size int/real/char", m_send%info(2:), " with comm ", &
          para_env%group, " and tag ", message_tag
      IF(m_send%info(2).GT.0) DEALLOCATE(m_send%task_int)
      IF(m_send%info(3).GT.0) DEALLOCATE(m_send%task_real)
      IF(m_send%info(4).GT.0) DEALLOCATE(m_send%task_char)
      IF(PRESENT(success)) success = .TRUE.
    END IF

    ! =============================
    ! broadcast
    ! =============================
    IF(dest.EQ.bcast_group) THEN
      IF(para_env%num_pe.GT.1) THEN
        CALL mp_bcast(m_send%info, MASTER_COMM_ID, para_env%group)
        IF(m_send%info(2).GT.0) THEN
          IF(.NOT.act_send_recv) ALLOCATE(m_send%task_int(m_send%info(2)))
          CALL mp_bcast(m_send%task_int, MASTER_COMM_ID, para_env%group)
        END IF
        IF(m_send%info(3).GT.0) THEN
          IF(.NOT.act_send_recv) ALLOCATE(m_send%task_real(m_send%info(3)))
          CALL mp_bcast(m_send%task_real, MASTER_COMM_ID, para_env%group)
        END IF
        IF(m_send%info(4).GT.0) THEN
          IF(.NOT.act_send_recv) ALLOCATE(m_send%task_char(m_send%info(3)))
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          !TODO bcast char CALL mp_bcast(m_send%task_char, MASTER_COMM_ID, para_env%group)
        END IF
      END IF
      ! sender delete arrays
      IF(act_send_recv) THEN
        IF(m_send%info(2).GT.0) DEALLOCATE(m_send%task_int)
        IF(m_send%info(3).GT.0) DEALLOCATE(m_send%task_real)
        IF(m_send%info(4).GT.0) DEALLOCATE(m_send%task_char)
      END IF
    END IF

    ! =============================
    ! receiving message
    ! =============================
    IF((act_send_recv.EQV.recv_msg).AND.dest.NE.bcast_group) THEN
      flag = .FALSE.
      tmp_tag = TMC_STATUS_WAIT_FOR_NEW_TASK
      IF(PRESENT(wait_for_message)) THEN
        dest = mp_any_source
        CALL mp_probe(dest, para_env%group, tmp_tag)
        flag = .TRUE.
      ELSE
        participant_loop: DO i=0, para_env%num_pe-1
          IF(i.NE.para_env%mepos) THEN
            dest = i
            CALL mp_probe(dest, para_env%group, tmp_tag)
            IF(dest.EQ.i)THEN
              flag = .TRUE.
              EXIT participant_loop
            END IF
          END IF
        END DO participant_loop
      END IF
      IF(flag.EQV..FALSE.) THEN
        IF(PRESENT(success)) success = .FALSE.
        DEALLOCATE(m_send)
        RETURN
      END IF

      IF(tmp_tag.EQ.TMC_STAT_SCF_STEP_ENER_RECEIVE) THEN
        ! CP2K send back SCF step energies without info message
        message_tag    = TMC_STAT_SCF_STEP_ENER_RECEIVE
        m_send%info(1) = TMC_STAT_SCF_STEP_ENER_RECEIVE
        m_send%info(2) = 0 ! no integer values
        m_send%info(3) = 1 ! one double values (SCF total energy)
        m_send%info(4) = 0 ! no character values
      ELSE
        message_tag = mp_any_tag
        ! first get message type and sizes
        CALL mp_recv(m_send%info, dest, message_tag, para_env%group)
      END IF
      IF(DEBUG.GE.1) &
        WRITE(tmc_out_file_nr,*)"TMC|message: ID: ",para_env%mepos, &
          " recv element info from ", dest, " of stat ", m_send%info(1), &
          " with size int/real/char", m_send%info(2:)
       !-- receive message integer part
      IF(m_send%info(2).GT.0) THEN
        ALLOCATE(m_send%task_int(m_send%info(2)))
        CALL mp_recv(m_send%task_int, dest, message_tag, para_env%group)
      END IF
      !-- receive message double (floatingpoint) part
      IF(m_send%info(3).GT.0) THEN
        ALLOCATE(m_send%task_real(m_send%info(3)))
        CALL mp_recv(m_send%task_real, dest, message_tag, para_env%group)
      END IF
      !-- receive message character part
      IF(m_send%info(4).GT.0) THEN
        ALLOCATE(m_send%task_char(m_send%info(4)))
        CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
        !TODO recv characters CALL mp_recv(m_send%task_char, dest, message_tag, para_env%group)
      END IF
    END IF

    ! handling recived message
    IF(act_send_recv.EQV.recv_msg) THEN
      ! if the element is supposed to be canceled but received message is not canceling receipt do not handle element
      ! (because element could be already deallocated, and hence a new element would be created -> not necessary)
      IF(PRESENT(elem_array))THEN
        IF(elem_array(dest)%canceled .AND. m_send%info(1).NE.TMC_CANCELING_RECEIPT) THEN
          msg_type = m_send%info(1)
          IF(m_send%info(2).GT.0) DEALLOCATE(m_send%task_int)
          IF(m_send%info(3).GT.0) DEALLOCATE(m_send%task_real)
          IF(m_send%info(4).GT.0) DEALLOCATE(m_send%task_char)
          ! to check for further messages
          IF(PRESENT(success)) success = .TRUE.
          DEALLOCATE(m_send)
          RETURN
        END IF
      END IF

      msg_type = m_send%info(1)
      SELECT CASE (m_send%info(1))
      CASE(TMC_STAT_START_CONF_REQUEST, TMC_CANCELING_MESSAGE, &
           TMC_CANCELING_RECEIPT, TMC_STATUS_WAIT_FOR_NEW_TASK, &
           TMC_STATUS_CALCULATING, TMC_STAT_ANALYSIS_RESULT)
        ! nothing to do here
      CASE(TMC_STATUS_WORKER_INIT)
        CALL read_worker_init_message(msg_type, tmc_params, m_send, error)
      CASE(TMC_STAT_START_CONF_RESULT, TMC_STAT_INIT_ANALYSIS)
        IF(PRESENT(elem_array)) THEN
          CALL read_start_conf_message(msg_type, elem_array(dest)%elem, &
                 result_count, m_send, tmc_params, error)
        ELSE
          CALL read_start_conf_message(msg_type, elem, result_count, m_send, &
                 tmc_params, error)
        END IF
      CASE(TMC_STAT_APPROX_ENERGY_RESULT)
        CALL read_approx_energy_result(elem_array(dest)%elem, m_send, tmc_params, error)
      CASE(TMC_STAT_ENERGY_REQUEST, TMC_STAT_APPROX_ENERGY_REQUEST)
        CALL read_energy_request_message(elem, m_send, tmc_params, error)
      CASE(TMC_STAT_ENERGY_RESULT)
        IF(PRESENT(elem_array)) &
          CALL read_energy_result_message(elem_array(dest)%elem, m_send, tmc_params, error)
      CASE(TMC_STAT_NMC_REQUEST, TMC_STAT_NMC_BROADCAST, &
           TMC_STAT_MD_REQUEST, TMC_STAT_MD_BROADCAST)
        CALL read_NMC_request_massage(msg_type, elem, m_send, tmc_params, error)
      CASE(TMC_STAT_NMC_RESULT, TMC_STAT_MD_RESULT)
        IF(PRESENT(elem_array)) &
          CALL read_NMC_result_massage(msg_type, elem_array(dest)%elem, m_send, tmc_params, error)
      CASE(TMC_STATUS_FAILED, TMC_STATUS_STOP_RECEIPT)
        ! if task is failed, handle situation in outer routine
      CASE(TMC_STAT_SCF_STEP_ENER_RECEIVE)
        CALL read_scf_step_ener(elem_array(dest)%elem, m_send, error)
      CASE(TMC_STAT_ANALYSIS_REQUEST)
        CALL read_analysis_request_message(elem, m_send, tmc_params, error)
      CASE DEFAULT
        CALL cp_assert(.FALSE., cp_failure_level,cp_assertion_failed,routineP,&
                       "try to receive unknown message type "//cp_to_string(msg_type)//&
                       "from source "//cp_to_string(dest),&
                       error, failure=flag)
      END SELECT
      IF(m_send%info(2).GT.0) DEALLOCATE(m_send%task_int)
      IF(m_send%info(3).GT.0) DEALLOCATE(m_send%task_real)
      IF(m_send%info(4).GT.0) DEALLOCATE(m_send%task_char)
      IF(PRESENT(success)) success = .TRUE.
    END IF

    ! ATTENTION there is also an short exit (RETURN) after probing for new messages
    DEALLOCATE(m_send)
  END SUBROUTINE tmc_message

! *****************************************************************************
!> \brief following pairs of routines creating outgoing and 
!>        handling incoming messages are defined
! *****************************************************************************

! *****************************************************************************
!> \brief set the messege just with an status tag
!> \param msg_type the status tag 
!> \param m_send the message structure
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE create_status_message(msg_type, m_send, error)
    INTEGER                                  :: msg_type
    TYPE(message_send), POINTER              :: m_send
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_status_message', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(m_send),cp_failure_level,routineP,error,failure)

    ! nothing to do, send just the message tag

    CPPostconditionNoFail(.NOT.ALLOCATED(m_send%task_int),cp_failure_level,routineP,error)
    CPPostconditionNoFail(.NOT.ALLOCATED(m_send%task_real),cp_failure_level,routineP,error)

  END SUBROUTINE create_status_message

  !============================================================================
  ! message for requesting start configuration
  !============================================================================
!! *****************************************************************************
!!> \brief the message for sending the atom mass 
!!>        (number of atoms is also tranfered)
!!>        atom names have to be done seperately, 
!!>        because character send only with bcast possible
!!> \param tmc_parms th send the cell properties
!!> \param m_send the message structure
!!> \param error variable to control error logging, stopping,...
!!>        see module cp_error_handling
!!> \author Mandes 02.2013
!! *****************************************************************************
!  SUBROUTINE create_atom_mass_message(m_send, atoms,error)
!    TYPE(tmc_atom_type), DIMENSION(:), POINTER    :: atoms
!    TYPE(message_send), POINTER              :: m_send
!    TYPE(cp_error_type), INTENT(inout)       :: error
!
!    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_atom_mass_message', &
!      routineP = moduleN//':'//routineN
!
!    INTEGER                                  :: counter, i, &
!                                                msg_size_real
!    LOGICAL                                  :: failure
!
!    failure = .FALSE.
!
!    CPPrecondition(ASSOCIATED(m_send),cp_failure_level,routineP,error,failure)
!    CPPrecondition(.NOT.ALLOCATED(m_send%task_int),cp_failure_level,routineP,error,failure)
!    CPPrecondition(.NOT.ALLOCATED(m_send%task_real),cp_failure_level,routineP,error,failure)
!    CPPrecondition(.NOT.ALLOCATED(m_send%task_char),cp_failure_level,routineP,error,failure)
!
!    counter =1
!    msg_size_real = 1+SIZE(tmc_params%cell%hmat)+ 1+SIZE(atoms) +1 
!    ALLOCATE(m_send%task_real(msg_size_real))
!
!    m_send%task_real(1) = REAL(SIZE(atoms,KIND=dp))
!    DO i=1, SIZE(atoms)
!      m_send%task_real(counter+i) = atoms(i)%mass
!    END DO
!    counter = counter + 1+INT(m_send%task_real(counter))
!    m_send%task_real(counter) = REAL(message_end_flag, KIND=dp) !message end
!    CPPostconditionNoFail(INT(m_send%task_real(msg_size_real)).EQ.message_end_flag,cp_failure_level,routineP,error)
!  END SUBROUTINE create_atom_mass_message
!
!! *****************************************************************************
!!> \brief the message for reading the atom mass 
!!>        (number of atoms is also tranfered)
!!>        atom names have to be done seperately, 
!!>        because character send only with bcast possible
!!> \param tmc_parms th send the cell properties
!!> \param m_send the message structure
!!> \param error variable to control error logging, stopping,...
!!>        see module cp_error_handling
!!> \author Mandes 02.2013
!! *****************************************************************************
!  SUBROUTINE read_atom_mass_message(m_send, atoms,error)
!    TYPE(tmc_atom_type), DIMENSION(:), &
!      POINTER                                :: atoms
!    TYPE(message_send), POINTER              :: m_send
!    TYPE(cp_error_type), INTENT(inout)       :: error
!
!    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_atom_mass_message', &
!      routineP = moduleN//':'//routineN
!
!    INTEGER                                  :: counter, i, nr_atoms
!    LOGICAL                                  :: failure
!
!    failure = .FALSE.
!
!    CPPrecondition(ASSOCIATED(m_send),cp_failure_level,routineP,error,failure)
!    CPPrecondition(.NOT.ALLOCATED(m_send%task_int),cp_failure_level,routineP,error,failure)
!    CPPrecondition(ALLOCATED(m_send%task_real),cp_failure_level,routineP,error,failure)
!    CPPrecondition(.NOT.ALLOCATED(m_send%task_char),cp_failure_level,routineP,error,failure)
!
!    counter =1
!    nr_atoms = m_send%task_real(counter)
!    IF(.NOT.ASSOCIATED(atoms)) CALL allocate_tmc_atom_type(atoms, nr_atoms, error)
!    DO i=1, SIZE(atoms)
!      atoms(i)%mass = m_send%task_real(counter+i)
!    END DO
!    counter = counter + 1+INT(m_send%task_real(counter))
!    CPPostconditionNoFail(INT(m_send%task_real(counter)).EQ.message_end_flag,cp_failure_level,routineP,error)
!  END SUBROUTINE read_atom_mass_message

! *****************************************************************************
!> \brief the message for the initial values (cell size) to the workers
!> \param msg_type the status tag 
!> \param tmc_parms to send the cell properties
!> \param m_send the message structure
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 07.2013
! *****************************************************************************
  SUBROUTINE create_worker_init_message(msg_type, tmc_params, m_send, error)
    INTEGER                                  :: msg_type
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(message_send), POINTER              :: m_send
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_worker_init_message', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter, msg_size_int, &
                                                msg_size_real
    LOGICAL                                  :: failure

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(m_send),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ALLOCATED(m_send%task_int),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ALLOCATED(m_send%task_real),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ALLOCATED(m_send%task_char),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params%cell),cp_failure_level,routineP,error,failure)
    
    counter = 1
    msg_size_int = 1+SIZE(tmc_params%cell%perd) +1+1 +1 +1
    ALLOCATE(m_send%task_int(msg_size_int))
    m_send%task_int(counter) = SIZE(tmc_params%cell%perd) ! periodicity of the cell
    counter = counter +1+m_send%task_int(counter)
    m_send%task_int(2:counter-1) = tmc_params%cell%perd(:)
    m_send%task_int(counter) = 1
    m_send%task_int(counter+1) = tmc_params%cell%symmetry_id
    m_send%task_int(counter+2) = 0
    IF(tmc_params%cell%orthorhombic) m_send%task_int(counter+2) = 1
    counter = counter +3
    m_send%task_int(counter) = message_end_flag
    CPPostconditionNoFail(counter.EQ.SIZE(m_send%task_int),cp_failure_level,routineP,error)

    !float array with cell vectors
    msg_size_real = 1+SIZE(tmc_params%cell%hmat) +1
    ALLOCATE(m_send%task_real(msg_size_real))
    counter = 1
    m_send%task_real(counter) = SIZE(tmc_params%cell%hmat) ! cell vectors for cell size
    m_send%task_real(counter+1:counter+SIZE(tmc_params%cell%hmat)) = &
                                    RESHAPE(tmc_params%cell%hmat(:,:),&
                                     (/SIZE(tmc_params%cell%hmat)/))
    counter = counter + 1+INT(m_send%task_real(counter))
    m_send%task_real(counter) = REAL(message_end_flag, KIND=dp) !message end
    CPPostconditionNoFail(SIZE(m_send%task_real).EQ.msg_size_real,cp_failure_level,routineP,error)
    CPPostconditionNoFail(INT(m_send%task_real(msg_size_real)).EQ.message_end_flag,cp_failure_level,routineP,error)
  END SUBROUTINE create_worker_init_message

! *****************************************************************************
!> \brief the message for the initial values (cell size) to the workers
!> \param msg_type the status tag 
!> \param tmc_parms to send the cell properties
!> \param m_send the message structure
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 07.2013
! *****************************************************************************
  SUBROUTINE read_worker_init_message(msg_type, tmc_params, m_send, error)
    INTEGER                                  :: msg_type
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(message_send), POINTER              :: m_send
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_worker_init_message', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter
    LOGICAL                                  :: failure, flag

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(m_send),cp_failure_level,routineP,error,failure)
    CPPrecondition(m_send%info(3).GE.4,cp_failure_level,routineP,error,failure)

    IF(.NOT.ASSOCIATED(tmc_params%cell)) ALLOCATE(tmc_params%cell)
    counter = 1
    !int array 
    flag = INT(m_send%task_int(1)) .EQ.SIZE(tmc_params%cell%perd)
    CPPrecondition(flag,cp_failure_level,routineP,error,failure)
    counter = 1+ m_send%task_int(1)+1
    tmc_params%cell%perd = m_send%task_int(2:counter-1)
    tmc_params%cell%symmetry_id = m_send%task_int(counter+1)
    tmc_params%cell%orthorhombic = .FALSE.
    IF(m_send%task_int(counter+2).EQ.1) tmc_params%cell%orthorhombic = .TRUE.
    counter = counter +3
    CPPostcondition(counter.EQ.m_send%info(2),cp_failure_level,routineP,error, failure)
    CPPostcondition(m_send%task_int(counter).EQ.message_end_flag,cp_failure_level,routineP,error, failure)

    !float array with cell vectors
    counter = 1
    flag = INT(m_send%task_real(counter)).EQ.SIZE(tmc_params%cell%hmat)
    CPPrecondition(flag,cp_failure_level,routineP,error,failure)
    tmc_params%cell%hmat = &
      RESHAPE(m_send%task_real(counter+1:counter+&
        SIZE(tmc_params%cell%hmat)),(/3,3/) )
    counter = counter + 1+INT(m_send%task_real(counter))

    CPPostcondition(counter.EQ.m_send%info(3),cp_failure_level,routineP,error, failure)
    CPPostconditionNoFail(INT(m_send%task_real(m_send%info(3))).EQ.message_end_flag,cp_failure_level,routineP,error)

  END SUBROUTINE read_worker_init_message

! *****************************************************************************
!> \brief the message for sending back the initial configuration
!> \param msg_type the status tag 
!> \param elem the intial tree element with initial coordinates and energy 
!>        (using the approximated potential)
!> \param tmc_parms to send the cell properties
!> \param m_send the message structure
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE create_start_conf_message(msg_type, elem, result_count, &
                                       tmc_params, m_send, error)
    INTEGER                                  :: msg_type
    TYPE(tree_type), POINTER                 :: elem
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: result_count
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(message_send), POINTER              :: m_send
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_start_conf_message', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter, i, msg_size_int, &
                                                msg_size_real
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(m_send),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params%atoms),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ALLOCATED(m_send%task_int),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ALLOCATED(m_send%task_real),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ALLOCATED(m_send%task_char),cp_failure_level,routineP,error,failure)

    counter = 1
    msg_size_int = 1+SIZE(tmc_params%cell%perd) +1+1 +1 +1+SIZE(elem%mol) +1
    IF(msg_type.EQ.TMC_STAT_INIT_ANALYSIS) THEN
      CPPrecondition(PRESENT(result_count),cp_failure_level,routineP,error,failure)
      CPPrecondition(ASSOCIATED(result_count),cp_failure_level,routineP,error,failure)
      msg_size_int = msg_size_int + 1+ SIZE(result_count(1:))
    END IF
    ALLOCATE(m_send%task_int(msg_size_int))
    m_send%task_int(counter) = SIZE(tmc_params%cell%perd) ! periodicity of the cell
    counter = counter +1+m_send%task_int(counter)
    m_send%task_int(2:counter-1) = tmc_params%cell%perd(:)
    m_send%task_int(counter) = 1
    m_send%task_int(counter+1) = tmc_params%cell%symmetry_id
    m_send%task_int(counter+2) = 0
    IF(tmc_params%cell%orthorhombic) m_send%task_int(counter+2) = 1
    counter = counter +3
    m_send%task_int(counter) = SIZE(elem%mol)
    m_send%task_int(counter+1:counter+m_send%task_int(counter)) = elem%mol(:)
    counter = counter +1+m_send%task_int(counter)
    IF(msg_type.EQ.TMC_STAT_INIT_ANALYSIS) THEN
      m_send%task_int(counter) = SIZE(result_count(1:))
      m_send%task_int(counter+1:counter+m_send%task_int(counter)) = &
        result_count(1:)
      counter = counter +1+m_send%task_int(counter)
    END IF
    m_send%task_int(counter) = message_end_flag
    CPPostconditionNoFail(counter.EQ.SIZE(m_send%task_int),cp_failure_level,routineP,error)

    counter = 0
    !float array with pos, cell vectors, atom_mass
    msg_size_real = 1+SIZE(elem%pos) +1+SIZE(tmc_params%cell%hmat) &
                   +1+SIZE(tmc_params%atoms) +1
    ALLOCATE(m_send%task_real(msg_size_real))
    m_send%task_real(1) = REAL(SIZE(elem%pos),KIND=dp) ! positions
    counter = 2 + m_send%task_real(1)
    m_send%task_real(2:counter-1) = elem%pos
    m_send%task_real(counter) = SIZE(tmc_params%cell%hmat) ! cell vectors for cell size
    m_send%task_real(counter+1:counter+SIZE(tmc_params%cell%hmat)) = &
                                    RESHAPE(tmc_params%cell%hmat(:,:),&
                                            (/SIZE(tmc_params%cell%hmat)/))
    counter = counter + 1+INT(m_send%task_real(counter))
    m_send%task_real(counter) = SIZE(tmc_params%atoms) ! atom mass
    DO i=1, SIZE(tmc_params%atoms)
      m_send%task_real(counter+i) = tmc_params%atoms(i)%mass
    END DO
    counter = counter + 1+INT(m_send%task_real(counter))
    m_send%task_real(counter) = REAL(message_end_flag, KIND=dp) !message end
    CPPostconditionNoFail(SIZE(m_send%task_real).EQ.msg_size_real,cp_failure_level,routineP,error)
    CPPostconditionNoFail(INT(m_send%task_real(msg_size_real)).EQ.message_end_flag,cp_failure_level,routineP,error)

  END SUBROUTINE create_start_conf_message

! *****************************************************************************
!> \brief the message for sending back the initial configuration
!> \param msg_type the status tag 
!> \param elem the intial tree element with initial coordinates and energy 
!>        (using the approximated potential)
!> \param m_send the message structure
!> \param tmc_params the param struct with necessary values for allocation
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE read_start_conf_message(msg_type, elem, result_count, m_send, &
                                     tmc_params, error)
    INTEGER                                  :: msg_type
    TYPE(tree_type), POINTER                 :: elem
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: result_count
    TYPE(message_send), POINTER              :: m_send
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_start_conf_message', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter, i
    LOGICAL                                  :: failure, flag

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(tmc_params%atoms),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(m_send),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(m_send%info(3).GE.4,cp_failure_level,routineP,error,failure)

    IF(.NOT.ASSOCIATED(tmc_params%cell)) ALLOCATE(tmc_params%cell)
    CALL allocate_new_sub_tree_node(tmc_params=tmc_params, next_el=elem, &
                                    nr_dim=NINT(m_send%task_real(1)), error=error)
    counter = 1
    !int array 
    flag = INT(m_send%task_int(1)) .EQ.SIZE(tmc_params%cell%perd)
    CPPrecondition(flag,cp_failure_level,routineP,error,failure)
    counter = 1+ m_send%task_int(1)+1
    tmc_params%cell%perd = m_send%task_int(2:counter-1)
    tmc_params%cell%symmetry_id = m_send%task_int(counter+1)
    tmc_params%cell%orthorhombic = .FALSE.
    IF(m_send%task_int(counter+2).EQ.1) tmc_params%cell%orthorhombic = .TRUE.
    counter = counter +3
    elem%mol(:) = m_send%task_int(counter+1:counter+m_send%task_int(counter))
    counter = counter +1+m_send%task_int(counter)
    IF(msg_type.EQ.TMC_STAT_INIT_ANALYSIS) THEN
      CPPrecondition(PRESENT(result_count),cp_failure_level,routineP,error,failure)
      CPPrecondition(.NOT.ASSOCIATED(result_count),cp_failure_level,routineP,error,failure)
      ALLOCATE(result_count(m_send%task_int(counter)))
      result_count(:) = m_send%task_int(counter+1:counter+m_send%task_int(counter))
      counter = counter +1+m_send%task_int(counter)
    END IF
    CPPostcondition(counter.EQ.m_send%info(2),cp_failure_level,routineP,error, failure)
    CPPostcondition(m_send%task_int(counter).EQ.message_end_flag,cp_failure_level,routineP,error, failure)

    counter = 0
    !float array with pos, cell vectors, atom_mass
    counter = 2 + INT(m_send%task_real(1))
    elem%pos = m_send%task_real(2:counter-1)
    flag = INT(m_send%task_real(counter)).EQ.SIZE(tmc_params%cell%hmat)
    CPPrecondition(flag,cp_failure_level,routineP,error,failure)
    tmc_params%cell%hmat = &
      RESHAPE(m_send%task_real(counter+1:counter+&
        SIZE(tmc_params%cell%hmat)),(/3,3/) )
    counter = counter + 1+INT(m_send%task_real(counter))

    CALL allocate_tmc_atom_type(atoms=tmc_params%atoms, &
                                nr_atoms=INT(m_send%task_real(counter)), &
                                error=error)
    DO i=1, SIZE(tmc_params%atoms)
      tmc_params%atoms(i)%mass = m_send%task_real(counter+i)
    END DO
    counter = counter + 1+INT(m_send%task_real(counter))

    CPPostcondition(counter.EQ.m_send%info(3),cp_failure_level,routineP,error, failure)
    CPPostconditionNoFail(INT(m_send%task_real(m_send%info(3))).EQ.message_end_flag,cp_failure_level,routineP,error)

  END SUBROUTINE read_start_conf_message

  !============================================================================
  ! Energy messages
  !============================================================================
! *****************************************************************************
!> \brief creating message for requesting exact energy of new configuration
!> \param msg_type the status tag 
!> \param elem tree element with new coordinates 
!> \param m_send the message structure
!> \param tmc_params stuct with parameters (global settings)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE create_energy_request_message(msg_type, elem, m_send, &
                                           tmc_params, error)
    INTEGER                                  :: msg_type
    TYPE(tree_type), POINTER                 :: elem
    TYPE(message_send), POINTER              :: m_send
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'create_energy_request_message', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter, msg_size_int, &
                                                msg_size_real
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(m_send),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ALLOCATED(m_send%task_int),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ALLOCATED(m_send%task_real),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

    counter = 0
    !first integer array 
    msg_size_int = 1+1 +1+1 +1 ! 1+SIZE(elem%sub_tree_nr) +1+SIZE(elem%nr)
    ALLOCATE(m_send%task_int(msg_size_int))
    counter = 1
    m_send%task_int(counter) = 1 !SIZE(elem%sub_tree_nr)
    m_send%task_int(counter+1:counter+m_send%task_int(counter)) = elem%sub_tree_nr
    counter = counter + 1+m_send%task_int(counter)
    m_send%task_int(counter) = 1 !SIZE(elem%nr)
    m_send%task_int(counter+1:counter+m_send%task_int(counter)) = elem%nr
    counter = counter + 1+m_send%task_int(counter)
    m_send%task_int(counter) = message_end_flag
    CPPostconditionNoFail(SIZE(m_send%task_int).EQ.msg_size_int,cp_failure_level,routineP,error)
    CPPostconditionNoFail(m_send%task_int(msg_size_int).EQ.message_end_flag,cp_failure_level,routineP,error)

    !then float array with pos
    msg_size_real = 1+SIZE(elem%pos) +1
    IF(tmc_params%pressure.GE.0.0_dp) msg_size_real = msg_size_real+ 1+SIZE(elem%box_scale(:))
    ALLOCATE(m_send%task_real(msg_size_real))
    m_send%task_real(1) = SIZE(elem%pos)
    counter = 2 + m_send%task_real(1)
    m_send%task_real(2:counter-1) = elem%pos
    IF(tmc_params%pressure.GE.0.0_dp) THEN
       m_send%task_real(counter) = SIZE(elem%box_scale)
       m_send%task_real(counter+1:counter+INT(m_send%task_real(counter))) = elem%box_scale(:)
       counter = counter + 1+m_send%task_real(counter)
    END IF
    m_send%task_real(counter) = REAL(message_end_flag, KIND=dp) !message end

    CPPostconditionNoFail(SIZE(m_send%task_real).EQ.msg_size_real,cp_failure_level,routineP,error)
    CPPostconditionNoFail(INT(m_send%task_real(msg_size_real)).EQ.message_end_flag,cp_failure_level,routineP,error)
  END SUBROUTINE create_energy_request_message

! *****************************************************************************
!> \brief reading message for requesting exact energy of new configuration
!> \param elem tree element with new coordinates 
!> \param m_send the message structure
!> \param tmc_params stuct with parameters (global settings)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE read_energy_request_message(elem, m_send, tmc_params, error)
    TYPE(tree_type), POINTER                 :: elem
    TYPE(message_send), POINTER              :: m_send
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_energy_request_message', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(m_send),cp_failure_level,routineP,error,failure)
    CPPrecondition(m_send%info(3).GT.0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(elem),cp_failure_level,routineP,error,failure)

    ! initialize the new sub tree element
    IF(.NOT.ASSOCIATED(elem))THEN
       CALL allocate_new_sub_tree_node(next_el=elem, nr_dim=NINT(m_send%task_real(1)),&
                                       tmc_params=tmc_params, error=error)
    END IF
    ! read the integer values
    CPPrecondition(m_send%info(2).GT.0,cp_failure_level,routineP,error,failure)
    counter = 1
    elem%sub_tree_nr = m_send%task_int(counter+1)
    counter = counter + 1+m_send%task_int(counter)
    elem%nr = m_send%task_int(counter+1)
    counter = counter + 1+m_send%task_int(counter)
    CPPostconditionNoFail(m_send%task_int(counter).EQ.message_end_flag,cp_failure_level,routineP,error)

    !float array with pos
    counter = 0
    counter = 1 + NINT(m_send%task_real(1))
    elem%pos = m_send%task_real(2:counter)
    counter = counter +1
    IF(tmc_params%pressure.GE.0.0_dp) THEN
       elem%box_scale(:) = m_send%task_real(counter+1:counter+INT(m_send%task_real(counter)))
       counter = counter +1+m_send%task_real(counter)
    END IF

    CPPostcondition(counter.EQ.m_send%info(3),cp_failure_level,routineP,error, failure)
    CPPostconditionNoFail(INT(m_send%task_real(m_send%info(3))).EQ.message_end_flag,cp_failure_level,routineP,error)
  END SUBROUTINE read_energy_request_message

! *****************************************************************************
!> \brief creating message for sending back the exact energy of new conf
!> \param msg_type the status tag
!> \param elem tree element  with calculated energy
!> \param m_send the message structure
!> \param tmc_params stuct with parameters (global settings)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE create_energy_result_message(msg_type, elem, m_send, tmc_params, error)
    INTEGER                                  :: msg_type
    TYPE(tree_type), POINTER                 :: elem
    TYPE(message_send), POINTER              :: m_send
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_energy_result_message', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter, msg_size_int, &
                                                msg_size_real
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(m_send),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ALLOCATED(m_send%task_int),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ALLOCATED(m_send%task_real),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

    counter = 0
    !first integer array
    msg_size_int = 0
    ! for checking the tree element mapping, send back the tree numbers
    IF(DEBUG.GT.0)THEN
       msg_size_int = 1+1 +1+1 +1 ! 1+SIZE(elem%sub_tree_nr) +1+SIZE(elem%nr)
       ALLOCATE(m_send%task_int(msg_size_int))
       counter = 1
       m_send%task_int(counter) = 1 !SIZE(elem%sub_tree_nr)
       m_send%task_int(counter+1:counter+m_send%task_int(counter)) = elem%sub_tree_nr
       counter = counter + 1+m_send%task_int(counter)
       m_send%task_int(counter) = 1 !SIZE(elem%nr)
       m_send%task_int(counter+1:counter+m_send%task_int(counter)) = elem%nr
       counter = counter + m_send%task_int(counter) +1
       m_send%task_int(counter) = message_end_flag !message end
    END IF

    !then float array with energy of exact potential
    msg_size_real = 1+1 +1
    IF(tmc_params%print_forces) msg_size_real = msg_size_real + 1+SIZE(elem%frc)
    IF(tmc_params%print_dipole) msg_size_real = msg_size_real + 1+SIZE(elem%dipole)

    ALLOCATE(m_send%task_real(msg_size_real))
    m_send%task_real(1) = 1
    m_send%task_real(2) = elem%potential
    counter = 3
    IF(tmc_params%print_forces) THEN
       m_send%task_real(counter) = SIZE(elem%frc)
       m_send%task_real(counter+1:counter+NINT(m_send%task_real(counter))) = elem%frc
       counter = counter + NINT(m_send%task_real(counter)) +1
    END IF
    IF(tmc_params%print_dipole)THEN
       m_send%task_real(counter) = SIZE(elem%dipole)
       m_send%task_real(counter+1:counter+NINT(m_send%task_real(counter))) = elem%dipole
       counter = counter + NINT(m_send%task_real(counter)) +1
    END IF

    m_send%task_real(counter) = REAL(message_end_flag, KIND=dp) !message end

    CPPostconditionNoFail(SIZE(m_send%task_real).EQ.msg_size_real,cp_failure_level,routineP,error)
    CPPostconditionNoFail(INT(m_send%task_real(msg_size_real)).EQ.message_end_flag,cp_failure_level,routineP,error)
  END SUBROUTINE create_energy_result_message

! *****************************************************************************
!> \brief reading message for sending back the exact energy of new conf
!> \param elem tree element for storing new energy
!> \param m_send the message structure
!> \param tmc_params stuct with parameters (global settings)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE read_energy_result_message(elem, m_send, tmc_params, error)
    TYPE(tree_type), POINTER                 :: elem
    TYPE(message_send), POINTER              :: m_send
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_energy_result_message', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(m_send),cp_failure_level,routineP,error,failure)
    CPPrecondition(m_send%info(3).GT.0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)
    
    ! read the integer values
    ! for checking the tree element mapping, check the tree numbers
    IF(DEBUG.GT.0)THEN 
       counter = 1
       IF(elem%sub_tree_nr .NE. m_send%task_int(counter+1) .OR. &
          elem%nr .NE. m_send%task_int(counter+3)) THEN
          WRITE(6,*)"ERROR: read_energy_result: master got energy result of subtree elem ", &
                    m_send%task_int(counter+1),  m_send%task_int(counter+3), &
                    " but expect result of subtree elem", elem%sub_tree_nr, elem%nr
          STOP "read_energy_result: got energy result from unexpected tree element."
       END IF
    ELSE
      CPPrecondition(m_send%info(2).EQ.0,cp_failure_level,routineP,error,failure)
    END IF

    !then float array with energy of exact potential
    elem%potential = m_send%task_real(2)
    counter = 3
    IF(tmc_params%print_forces) THEN
       elem%frc(:) = m_send%task_real((counter+1):(counter+NINT(m_send%task_real(counter))))
       counter = counter + 1+NINT(m_send%task_real(counter))   
    END IF
    IF(tmc_params%print_dipole)THEN
      elem%dipole(:) = m_send%task_real((counter+1):(counter+NINT(m_send%task_real(counter))))
      counter = counter + 1+NINT(m_send%task_real(counter))
    END IF

    CPPostcondition(counter.EQ.m_send%info(3),cp_failure_level,routineP,error, failure)
    CPPostconditionNoFail(INT(m_send%task_real(m_send%info(3))).EQ.message_end_flag,cp_failure_level,routineP,error)
  END SUBROUTINE read_energy_result_message

! *****************************************************************************
!> \brief create message for sending back the approximate energy of new conf
!> \param elem tree element with calculated approx energy 
!> \param m_send the message structure
!> \param tmc_params stuct with parameters (global settings)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE create_approx_energy_result_message(msg_type, elem, m_send, &
                                                 tmc_params, error)
    INTEGER                                  :: msg_type
    TYPE(tree_type), POINTER                 :: elem
    TYPE(message_send), POINTER              :: m_send
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'create_approx_energy_result_message', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter, msg_size_real
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(m_send),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ALLOCATED(m_send%task_int),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ALLOCATED(m_send%task_real),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

    counter = 0

    !then float array with energy of exact potential
    msg_size_real = 1+1 +1
    IF(tmc_params%pressure.GE.0.0_dp) msg_size_real = msg_size_real+ 1+SIZE(elem%box_scale(:))

    ALLOCATE(m_send%task_real(msg_size_real))
    m_send%task_real(1) = 1
    m_send%task_real(2) = elem%e_pot_approx
    counter = 3
    ! the box size for NpT
    IF(tmc_params%pressure.GE.0.0_dp) THEN
       m_send%task_real(counter) = SIZE(elem%box_scale)
       m_send%task_real(counter+1:counter+INT(m_send%task_real(counter))) = elem%box_scale(:)
       counter = counter + 1+m_send%task_real(counter)
    END IF
    m_send%task_real(counter) = REAL(message_end_flag, KIND=dp) !message end

    CPPostconditionNoFail(SIZE(m_send%task_real).EQ.msg_size_real,cp_failure_level,routineP,error)
    CPPostconditionNoFail(INT(m_send%task_real(msg_size_real)).EQ.message_end_flag,cp_failure_level,routineP,error)
  END SUBROUTINE create_approx_energy_result_message

! *****************************************************************************
!> \brief reading message for sending back the exact energy of new conf
!> \param elem tree element for storing new energy
!> \param m_send the message structure
!> \param tmc_params the param struct with necessary parameters
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE read_approx_energy_result(elem, m_send, tmc_params, error)
    TYPE(tree_type), POINTER                 :: elem
    TYPE(message_send), POINTER              :: m_send
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_approx_energy_result', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(m_send),cp_failure_level,routineP,error,failure)
    CPPrecondition(m_send%info(2).EQ.0.AND.m_send%info(3).GT.0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

    !then float array with energy of exact potential
    elem%e_pot_approx = m_send%task_real(2)
    counter = 3
    IF(tmc_params%pressure.GE.0.0_dp) THEN
       elem%box_scale(:) = m_send%task_real(counter+1:counter+INT(m_send%task_real(counter)))
       counter = counter +1+m_send%task_real(counter)
    END IF

    CPPostcondition(counter.EQ.m_send%info(3),cp_failure_level,routineP,error, failure)
    CPPostconditionNoFail(INT(m_send%task_real(m_send%info(3))).EQ.message_end_flag,cp_failure_level,routineP,error)
  END SUBROUTINE read_approx_energy_result

  !============================================================================
  ! Nested Monte Carlo request messages
  !============================================================================
! *****************************************************************************
!> \brief creating message for Nested Monte Carlo sampling of new configuration
!> \param msg_type the status tag
!> \param elem tree element  with calculated energy
!> \param m_send the message structure
!> \param tmc_params stuct with parameters (global settings)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE create_NMC_request_massage(msg_type, elem, m_send, &
                                        tmc_params, error)
    INTEGER                                  :: msg_type
    TYPE(tree_type), POINTER                 :: elem
    TYPE(message_send), POINTER              :: m_send
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_NMC_request_massage', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter, msg_size_int, &
                                                msg_size_real
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(m_send),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ALLOCATED(m_send%task_int),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ALLOCATED(m_send%task_real),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

    counter = 0
    !first integer array with element status,mol_info, move type, sub tree, element nr, temp index
    msg_size_int = 1+SIZE(elem%elem_stat) +1+SIZE(elem%mol) +1+1 +1+1 +1+1 +1+1 +1

    ALLOCATE(m_send%task_int(msg_size_int))
    ! element status
    m_send%task_int(1) = SIZE(elem%elem_stat)
    counter = 2 + m_send%task_int(1)
    m_send%task_int(2:counter-1) = elem%elem_stat
    m_send%task_int(counter) = SIZE(elem%mol)
    m_send%task_int(counter+1:counter+m_send%task_int(counter)) = elem%mol(:)
    counter = counter +1+m_send%task_int(counter)
    ! element move type
    m_send%task_int(counter) = 1
    m_send%task_int(counter+1) = elem%move_type
    counter = counter +2 
    m_send%task_int(counter) = 1
    m_send%task_int(counter+1) = elem%nr
    counter = counter +2
    m_send%task_int(counter) = 1
    m_send%task_int(counter+1) = elem%sub_tree_nr
    counter = counter +2
    m_send%task_int(counter) = 1
    m_send%task_int(counter+1) = elem%temp_created
    m_send%task_int(counter+2) = message_end_flag !message end

    counter = 0
    !then float array with pos, (vel), random number seed, subbox_center
    msg_size_real = 1+SIZE(elem%pos) +1+SIZE(elem%rng_seed) +1+SIZE(elem%subbox_center(:)) +1
    IF(msg_type.EQ.TMC_STAT_MD_REQUEST.OR.msg_type.EQ.TMC_STAT_MD_BROADCAST) &
       msg_size_real = msg_size_real +1+SIZE(elem%vel)                                ! the velocities
    IF(tmc_params%pressure.GE.0.0_dp) msg_size_real = msg_size_real+ 1+SIZE(elem%box_scale(:)) ! box size for NpT

    ALLOCATE(m_send%task_real(msg_size_real))
    m_send%task_real(1) = SIZE(elem%pos)
    counter = 2 + m_send%task_real(1)
    m_send%task_real(2:counter-1) = elem%pos
    IF(msg_type.EQ.TMC_STAT_MD_REQUEST.OR.msg_type.EQ.TMC_STAT_MD_BROADCAST) THEN
       m_send%task_real(counter) = SIZE(elem%vel)
       m_send%task_real(counter+1:counter+NINT(m_send%task_real(counter))) = elem%vel
       counter = counter + 1+NINT(m_send%task_real(counter))
    END IF
    ! rng seed
    m_send%task_real(counter) = SIZE(elem%rng_seed)
    m_send%task_real(counter+1:counter+SIZE(elem%rng_seed)) = RESHAPE(elem%rng_seed(:,:,:),(/SIZE(elem%rng_seed)/))
    counter = counter + NINT(m_send%task_real(counter)) +1
    ! sub box center
    m_send%task_real(counter) = SIZE(elem%subbox_center(:))
    m_send%task_real(counter+1:counter+SIZE(elem%subbox_center)) = elem%subbox_center(:)
    counter = counter +1+NINT(m_send%task_real(counter))
    ! the box size for NpT
    IF(tmc_params%pressure.GE.0.0_dp) THEN
       m_send%task_real(counter) = SIZE(elem%box_scale)
       m_send%task_real(counter+1:counter+INT(m_send%task_real(counter))) = elem%box_scale(:)
       counter = counter + 1+m_send%task_real(counter)
    END IF
    m_send%task_real(counter) = message_end_flag !message end

    CPPostconditionNoFail(SIZE(m_send%task_int).EQ.msg_size_int,cp_failure_level,routineP,error)
    CPPostconditionNoFail(SIZE(m_send%task_real).EQ.msg_size_real,cp_failure_level,routineP,error)
    CPPostconditionNoFail(m_send%task_int(msg_size_int).EQ.message_end_flag,cp_failure_level,routineP,error)
    CPPostconditionNoFail(INT(m_send%task_real(msg_size_real)).EQ.message_end_flag,cp_failure_level,routineP,error)
  END SUBROUTINE create_NMC_request_massage

! *****************************************************************************
!> \brief reading message for Nested Monte Carlo sampling of new configuration
!> \param msg_type the status tag
!> \param elem tree element with new coordinates 
!> \param m_send the message structure
!> \param tmc_params stuct with parameters (global settings)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE read_NMC_request_massage(msg_type, elem, m_send, &
                                      tmc_params, error)
    INTEGER                                  :: msg_type
    TYPE(tree_type), POINTER                 :: elem
    TYPE(message_send), POINTER              :: m_send
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_NMC_request_massage', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter, num_dim, &
                                                rnd_seed_size
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(m_send),cp_failure_level,routineP,error,failure)
    CPPrecondition(m_send%info(2).GT.5.AND.m_send%info(3).GT.8,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

    counter = 0
    !first integer array with number of dimentions and random seed size
    rnd_seed_size = m_send%task_int(1+m_send%task_int(1)+1)
 
    IF(.NOT.ASSOCIATED(elem))THEN
      CALL allocate_new_sub_tree_node(next_el=elem, nr_dim=NINT(m_send%task_real(1)), &
                                      tmc_params=tmc_params, error=error)
    END IF
    ! element status
    counter = 2 + m_send%task_int(1)
    elem%elem_stat = m_send%task_int(2:counter-1)
    elem%mol(:) = m_send%task_int(counter+1:counter+m_send%task_int(counter))
    counter = counter +1+m_send%task_int(counter)
    ! element move type
    elem%move_type = m_send%task_int(counter+1)
    counter = counter +2
    elem%nr = m_send%task_int(counter+1)
    counter = counter +2
    elem%sub_tree_nr = m_send%task_int(counter+1)
    counter = counter +2
    elem%temp_created = m_send%task_int(counter+1)
    counter = counter +2
    CPPostcondition(counter.EQ.m_send%info(2),cp_failure_level,routineP,error, failure)
 
    counter = 0
    !then float array with pos, (vel), subbox_center and temp
    num_dim = NINT(m_send%task_real(1))
    counter = 2 + m_send%task_real(1)
    elem%pos = m_send%task_real(2:counter-1)
    IF(msg_type.EQ.TMC_STAT_MD_REQUEST.OR.msg_type.EQ.TMC_STAT_MD_BROADCAST)THEN
       elem%vel = m_send%task_real(counter+1:counter+NINT(m_send%task_real(counter)))
       counter = counter + NINT(m_send%task_real(counter))+1
    END IF
    ! rng seed
    elem%rng_seed(:,:,:) = RESHAPE(m_send%task_real(counter+1:counter+SIZE(elem%rng_seed)), (/3,2,3/))
    counter = counter + NINT(m_send%task_real(counter))+1
    ! sub box center
    elem%subbox_center(:) = m_send%task_real(counter+1:counter+INT(m_send%task_real(counter))) 
    counter = counter +1+NINT(m_send%task_real(counter))

    IF(tmc_params%pressure .GE. 0.0_dp) THEN
       elem%box_scale(:) = m_send%task_real(counter+1:counter+INT(m_send%task_real(counter)))
       counter = counter +1+m_send%task_real(counter)
    ELSE 
       elem%box_scale(:) = 1.0_dp
    END IF

    CPPostcondition(counter.EQ.m_send%info(3),cp_failure_level,routineP,error, failure)
    CPPostconditionNoFail(m_send%task_int(m_send%info(2)).EQ.message_end_flag,cp_failure_level,routineP,error)
    CPPostconditionNoFail(INT(m_send%task_real(m_send%info(3))).EQ.message_end_flag,cp_failure_level,routineP,error)
  END SUBROUTINE read_NMC_request_massage

  !============================================================================
  ! Nested Monte Carlo RESULT messages
  !============================================================================
! *****************************************************************************
!> \brief creating message for Nested Monte Carlo sampling result
!> \param msg_type the status tag
!> \param elem tree element  with calculated energy
!> \param m_send the message structure
!> \param tmc_params environment with move types and sizes
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE create_NMC_result_massage(msg_type, elem, m_send, tmc_params,&
                                       error)
    INTEGER                                  :: msg_type
    TYPE(tree_type), POINTER                 :: elem
    TYPE(message_send), POINTER              :: m_send
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_NMC_result_massage', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter, msg_size_int, &
                                                msg_size_real
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(m_send),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ALLOCATED(m_send%task_int),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ALLOCATED(m_send%task_real),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

    !first integer array with status, nmc_acc_counts, subbox_acc_count and (subbox rejectance)
    msg_size_int = 1+SIZE(elem%mol) &
                  +1+SIZE(tmc_params%nmc_move_types%mv_count) &
                  +1+SIZE(tmc_params%nmc_move_types%acc_count) +1
    IF(DEBUG.GT.0) msg_size_int = msg_size_int + 1+1 +1+1
    IF(.NOT.ANY(tmc_params%sub_box_size.LE.0.1_dp)) &
       msg_size_int = msg_size_int +1+SIZE(tmc_params%nmc_move_types%subbox_count) &
                      +1+SIZE(tmc_params%nmc_move_types%subbox_acc_count)

    ALLOCATE(m_send%task_int(msg_size_int))
    counter = 1
    IF(DEBUG.GT.0) THEN
       ! send the element number back
       m_send%task_int(counter) = 1
       m_send%task_int(counter+1) = elem%sub_tree_nr
       counter = counter +1 +m_send%task_int(counter)
       m_send%task_int(counter) = 1
       m_send%task_int(counter+1) = elem%nr
       counter = counter +1 +m_send%task_int(counter)
    END IF
    ! the molecule information
    m_send%task_int(counter) = SIZE(elem%mol)
    m_send%task_int(counter+1:counter+m_send%task_int(counter)) = elem%mol(:)
    counter = counter +1+m_send%task_int(counter)
    ! the counters for each move type 
    m_send%task_int(counter) = SIZE(tmc_params%nmc_move_types%mv_count)
    m_send%task_int(counter+1:counter+m_send%task_int(counter)) = &
       RESHAPE(tmc_params%nmc_move_types%mv_count(:,:),&
               (/SIZE(tmc_params%nmc_move_types%mv_count)/))
    counter = counter +1 +m_send%task_int(counter)
    ! the counter for the accepted moves
    m_send%task_int(counter) = SIZE(tmc_params%nmc_move_types%acc_count)
    m_send%task_int(counter+1:counter+m_send%task_int(counter)) = &
       RESHAPE(tmc_params%nmc_move_types%acc_count(:,:),&
               (/SIZE(tmc_params%nmc_move_types%acc_count)/))
    counter = counter +1 +m_send%task_int(counter)
    ! amount of rejected subbox moves
    IF(.NOT.ANY(tmc_params%sub_box_size.LE.0.1_dp)) THEN
       m_send%task_int(counter) = SIZE(tmc_params%nmc_move_types%subbox_count)
       m_send%task_int(counter+1:counter+m_send%task_int(counter)) = &
          RESHAPE(tmc_params%nmc_move_types%subbox_count(:,:),&
                  (/SIZE(tmc_params%nmc_move_types%subbox_count)/))
       counter = counter +1 +m_send%task_int(counter)
       m_send%task_int(counter) = SIZE(tmc_params%nmc_move_types%subbox_acc_count)
       m_send%task_int(counter+1:counter+m_send%task_int(counter)) = &
          RESHAPE(tmc_params%nmc_move_types%subbox_acc_count(:,:),&
                  (/SIZE(tmc_params%nmc_move_types%subbox_acc_count)/))
       counter = counter +1 +m_send%task_int(counter)
    END IF
    m_send%task_int(counter) = message_end_flag ! message end

    counter = 0
    !then float array with pos,(vel, e_kin_befor_md, ekin),(forces),rng_seed,
    !                       potential,e_pot_approx,acc_prob,subbox_prob
    msg_size_real =  1+SIZE(elem%pos) &             ! pos
                    +1+SIZE(elem%rng_seed)&         ! rng_seed
                    +1+1 &                          ! potential
                    +1+1 &                          ! e_pot_approx
                    +1                              ! check bit

    IF(msg_type.EQ.TMC_STAT_MD_REQUEST.OR.msg_type.EQ.TMC_STAT_MD_RESULT.OR.&
       msg_type.EQ.TMC_STAT_MD_BROADCAST) &
       msg_size_real = msg_size_real +1+SIZE(elem%vel) +1+1 +1+1 ! for MD also: vel, e_kin_befor_md, ekin

    ALLOCATE(m_send%task_real(msg_size_real))
    ! pos
    counter = 1
    m_send%task_real(counter) = SIZE(elem%pos)
    m_send%task_real(counter+1:counter+NINT(m_send%task_real(counter))) = elem%pos
    counter = counter +1 +NINT(m_send%task_real(counter))
    ! rng seed
    m_send%task_real(counter) = SIZE(elem%rng_seed)
    m_send%task_real(counter+1:counter+SIZE(elem%rng_seed)) = &
                     RESHAPE(elem%rng_seed(:,:,:),(/SIZE(elem%rng_seed)/))
    counter = counter +1 +NINT(m_send%task_real(counter))
    ! potential
    m_send%task_real(counter) = 1
    m_send%task_real(counter+1) = elem%potential
    counter = counter +2
    ! approximate potential energy
    m_send%task_real(counter) = 1
    m_send%task_real(counter+1) = elem%e_pot_approx
    counter = counter +2
    ! for MD also: vel, e_kin_befor_md, ekin
    IF(msg_type.EQ.TMC_STAT_MD_REQUEST.OR.msg_type.EQ.TMC_STAT_MD_RESULT.OR.&
       msg_type.EQ.TMC_STAT_MD_BROADCAST)THEN
       m_send%task_real(counter) = SIZE(elem%vel)
       m_send%task_real(counter+1:counter+NINT(m_send%task_real(counter))) = elem%vel
       counter = counter + 1+m_send%task_real(counter)
       m_send%task_real(counter) = 1
       m_send%task_real(counter+1) = elem%ekin_before_md
       counter = counter+2
       m_send%task_real(counter) = 1
       m_send%task_real(counter+1) = elem%ekin
       counter = counter+2
    END IF
    m_send%task_real(counter) = message_end_flag ! message end

    CPPostconditionNoFail(SIZE(m_send%task_int).EQ.msg_size_int,cp_failure_level,routineP,error)
    CPPostconditionNoFail(SIZE(m_send%task_real).EQ.msg_size_real,cp_failure_level,routineP,error)
    CPPostconditionNoFail(m_send%task_int(msg_size_int).EQ.message_end_flag,cp_failure_level,routineP,error)
    CPPostconditionNoFail(INT(m_send%task_real(msg_size_real)).EQ.message_end_flag,cp_failure_level,routineP,error)
  END SUBROUTINE create_NMC_result_massage

! *****************************************************************************
!> \brief reading message for Nested Monte Carlo sampling result
!> \param msg_type the status tag
!> \param elem tree element  with calculated energy
!> \param m_send the message structure
!> \param tmc_params environment with move types and sizes
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE read_NMC_result_massage(msg_type, elem, m_send, tmc_params, &
                                     error)
    INTEGER                                  :: msg_type
    TYPE(tree_type), POINTER                 :: elem
    TYPE(message_send), POINTER              :: m_send
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_NMC_result_massage', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter
    INTEGER, DIMENSION(:, :), POINTER        :: acc_counter, mv_counter, &
                                                subbox_acc_counter, &
                                                subbox_counter
    LOGICAL                                  :: failure

    failure = .FALSE.
    NULLIFY(mv_counter, subbox_counter, acc_counter, subbox_acc_counter)

    CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(m_send),cp_failure_level,routineP,error,failure)
    CPPrecondition(m_send%info(2).GT.0.AND.m_send%info(3).GT.0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

    !first integer array with element status, random number seed, and move type
    counter = 1
    IF(DEBUG.GT.0) THEN
       IF((m_send%task_int(counter+1) .NE. elem%sub_tree_nr) .AND. (m_send%task_int(counter+3).NE.elem%nr)) &
          STOP "read_NMC_result_massage: got result of wrong element"
       counter = counter + 2+2
    END IF
    ! the molecule information
    elem%mol(:) = m_send%task_int(counter+1:counter+m_send%task_int(counter))
    counter = counter +1+m_send%task_int(counter)
    ! the counters for each move type
    ALLOCATE(mv_counter(0:SIZE(tmc_params%nmc_move_types%mv_count(:,1))-1,&
             SIZE(tmc_params%nmc_move_types%mv_count(1,:))))
    mv_counter(:,:) = RESHAPE(m_send%task_int(counter+1:counter+m_send%task_int(counter)), &
       (/SIZE(tmc_params%nmc_move_types%mv_count(:,1)), &
        SIZE(tmc_params%nmc_move_types%mv_count(1,:))/))
    counter = counter +1 +m_send%task_int(counter)
    ! the counter for the accepted moves
    ALLOCATE(acc_counter(0:SIZE(tmc_params%nmc_move_types%acc_count(:,1))-1,&
             SIZE(tmc_params%nmc_move_types%acc_count(1,:))))
    acc_counter(:,:) = RESHAPE(m_send%task_int(counter+1:counter+m_send%task_int(counter)), &
       (/SIZE(tmc_params%nmc_move_types%acc_count(:,1)), &
         SIZE(tmc_params%nmc_move_types%acc_count(1,:))/))
    counter = counter +1 +m_send%task_int(counter)
    ! amount of rejected subbox moves
    IF(.NOT.ANY(tmc_params%sub_box_size.LE.0.1_dp)) THEN
       ALLOCATE(subbox_counter(SIZE(tmc_params%nmc_move_types%subbox_count(:,1)), &
                SIZE(tmc_params%nmc_move_types%subbox_count(1,:))))
       subbox_counter(:,:) = RESHAPE(m_send%task_int(counter+1:counter+m_send%task_int(counter)), &
          (/SIZE(tmc_params%nmc_move_types%subbox_count(:,1)), &
            SIZE(tmc_params%nmc_move_types%subbox_count(1,:))/))
       counter = counter +1 +m_send%task_int(counter)
       ALLOCATE(subbox_acc_counter(SIZE(tmc_params%nmc_move_types%subbox_acc_count(:,1)), &
                SIZE(tmc_params%nmc_move_types%subbox_acc_count(1,:))))
       subbox_acc_counter(:,:) = RESHAPE(m_send%task_int(counter+1:counter+m_send%task_int(counter)), &
          (/SIZE(tmc_params%nmc_move_types%subbox_acc_count(:,1)), &
            SIZE(tmc_params%nmc_move_types%subbox_acc_count(1,:))/))
       counter = counter +1 +m_send%task_int(counter)   
    END IF
    CPPostcondition(counter.EQ.m_send%info(2),cp_failure_level,routineP,error, failure)
    
    counter = 0
    !then float array with pos, (vel, e_kin_befor_md, ekin), (forces), rng_seed, potential, e_pot_approx 
    counter = 1
    ! pos
    elem%pos = m_send%task_real(counter+1:counter+NINT(m_send%task_real(counter)))
    counter = counter +1 +NINT(m_send%task_real(counter))
    ! rng seed
    elem%rng_seed(:,:,:) = RESHAPE(m_send%task_real(counter+1:counter+SIZE(elem%rng_seed)), (/3,2,3/))
    counter = counter +1 +NINT(m_send%task_real(counter))
    ! potential
    elem%potential = m_send%task_real(counter+1)
    counter = counter +2
    ! approximate potential energy
    elem%e_pot_approx  = m_send%task_real(counter+1)
    counter = counter +2
    ! for MD also: vel, e_kin_befor_md, ekin
    IF(msg_type.EQ.TMC_STAT_MD_REQUEST.OR.msg_type.EQ.TMC_STAT_MD_RESULT.OR.&
       msg_type.EQ.TMC_STAT_MD_BROADCAST)THEN
       elem%vel = m_send%task_real(counter+1:counter+NINT(m_send%task_real(counter))) 
       counter = counter + 1+m_send%task_real(counter)
       IF(.NOT.(tmc_params%task_type.EQ.task_type_gaussian_adaptation)) &
          elem%ekin_before_md = m_send%task_real(counter+1)
       counter = counter+2
       elem%ekin = m_send%task_real(counter+1)
       counter = counter+2
    END IF

    CALL add_mv_prob(move_types=tmc_params%nmc_move_types, prob_opt=tmc_params%esimate_acc_prob, &
                     mv_counter=mv_counter,acc_counter=acc_counter, error=error)
    IF(.NOT.ANY(tmc_params%sub_box_size.LE.0.1_dp)) THEN
      CALL add_mv_prob(move_types=tmc_params%nmc_move_types, prob_opt=tmc_params%esimate_acc_prob, &
                       subbox_counter=subbox_counter, subbox_acc_counter=subbox_acc_counter, error=error)
    END IF

    DEALLOCATE(mv_counter, acc_counter)
    IF(.NOT.ANY(tmc_params%sub_box_size.LE.0.1_dp)) &
      DEALLOCATE(subbox_counter, subbox_acc_counter)
    CPPostcondition(counter.EQ.m_send%info(3),cp_failure_level,routineP,error, failure)
    CPPostconditionNoFail(m_send%task_int(m_send%info(2)).EQ.message_end_flag,cp_failure_level,routineP,error)
    CPPostconditionNoFail(INT(m_send%task_real(m_send%info(3))).EQ.message_end_flag,cp_failure_level,routineP,error)
  END SUBROUTINE read_NMC_result_massage

  !============================================================================
  ! Analysis element messages
  !============================================================================
! *****************************************************************************
!> \brief creating message for requesting analysing a new configuration
!>        we plot temperatur index into the sub tree number and 
!>        the Markov chain number into the element number
!> \param elem tree element with new coordinates 
!> \param m_send the message structure
!> \param tmc_params stuct with parameters (global settings)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE create_analysis_request_message(list_elem, m_send, &
                                           tmc_params, error)
    TYPE(elem_list_type), POINTER            :: list_elem
    TYPE(message_send), POINTER              :: m_send
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'create_analysis_request_message', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter, msg_size_int, &
                                                msg_size_real
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(m_send),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ALLOCATED(m_send%task_int),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ALLOCATED(m_send%task_real),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(list_elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

    counter = 0
    !first integer array 
    msg_size_int = 1+1 +1+1 +1 ! 1+SIZE(list_elem%elem%sub_tree_nr) +1+SIZE(list_elem%elem%nr)
    ALLOCATE(m_send%task_int(msg_size_int))
    counter = 1
    m_send%task_int(counter) = 1 ! temperature index 
    m_send%task_int(counter+1:counter+m_send%task_int(counter)) = list_elem%temp_ind
    counter = counter + 1+m_send%task_int(counter)
    m_send%task_int(counter) = 1 ! Markov chain number
    m_send%task_int(counter+1:counter+m_send%task_int(counter)) = list_elem%nr
    counter = counter + 1+m_send%task_int(counter)
    m_send%task_int(counter) = message_end_flag
    CPPostconditionNoFail(SIZE(m_send%task_int).EQ.msg_size_int,cp_failure_level,routineP,error)
    CPPostconditionNoFail(m_send%task_int(msg_size_int).EQ.message_end_flag,cp_failure_level,routineP,error)

    !then float array with pos
    msg_size_real = 1+SIZE(list_elem%elem%pos) +1
    IF(tmc_params%pressure.GE.0.0_dp) msg_size_real = msg_size_real+ 1+SIZE(list_elem%elem%box_scale(:))
    ALLOCATE(m_send%task_real(msg_size_real))
    m_send%task_real(1) = SIZE(list_elem%elem%pos)
    counter = 2 + m_send%task_real(1)
    m_send%task_real(2:counter-1) = list_elem%elem%pos
    IF(tmc_params%pressure.GE.0.0_dp) THEN
       m_send%task_real(counter) = SIZE(list_elem%elem%box_scale)
       m_send%task_real(counter+1:counter+INT(m_send%task_real(counter))) = list_elem%elem%box_scale(:)
       counter = counter + 1+m_send%task_real(counter)
    END IF
    m_send%task_real(counter) = REAL(message_end_flag, KIND=dp) !message end

    CPPostconditionNoFail(SIZE(m_send%task_real).EQ.msg_size_real,cp_failure_level,routineP,error)
    CPPostconditionNoFail(INT(m_send%task_real(msg_size_real)).EQ.message_end_flag,cp_failure_level,routineP,error)
  END SUBROUTINE create_analysis_request_message

! *****************************************************************************
!> \brief reading message for requesting exact energy of new configuration
!> \param elem tree element with new coordinates 
!> \param m_send the message structure
!> \param tmc_params stuct with parameters (global settings)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE read_analysis_request_message(elem, m_send, tmc_params, error)
    TYPE(tree_type), POINTER                 :: elem
    TYPE(message_send), POINTER              :: m_send
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'read_analysis_request_message', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(m_send),cp_failure_level,routineP,error,failure)
    CPPrecondition(m_send%info(3).GT.0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(elem),cp_failure_level,routineP,error,failure)

    ! initialize the new sub tree element
    IF(.NOT.ASSOCIATED(elem))THEN
       CALL allocate_new_sub_tree_node(next_el=elem, nr_dim=NINT(m_send%task_real(1)),&
                                       tmc_params=tmc_params, error=error)
    END IF
    ! read the integer values
    CPPrecondition(m_send%info(2).GT.0,cp_failure_level,routineP,error,failure)
    counter = 1
    elem%sub_tree_nr = m_send%task_int(counter+1)
    counter = counter + 1+m_send%task_int(counter)
    elem%nr = m_send%task_int(counter+1)
    counter = counter + 1+m_send%task_int(counter)
    CPPostconditionNoFail(m_send%task_int(counter).EQ.message_end_flag,cp_failure_level,routineP,error)

    !float array with pos
    counter = 0
    counter = 1 + NINT(m_send%task_real(1))
    elem%pos = m_send%task_real(2:counter)
    counter = counter +1
    IF(tmc_params%pressure.GE.0.0_dp) THEN
       elem%box_scale(:) = m_send%task_real(counter+1:counter+INT(m_send%task_real(counter)))
       counter = counter +1+m_send%task_real(counter)
    END IF

    CPPostcondition(counter.EQ.m_send%info(3),cp_failure_level,routineP,error, failure)
    CPPostconditionNoFail(INT(m_send%task_real(m_send%info(3))).EQ.message_end_flag,cp_failure_level,routineP,error)
  END SUBROUTINE read_analysis_request_message

  !============================================================================
  ! SCF step energies (receiving from CP2K)
  !============================================================================
! *****************************************************************************
!> \brief routine cancel the other group participants
!> \param elem tree element  with approximated energy
!> \param m_send the message structure
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE read_scf_step_ener(elem, m_send, error)
    TYPE(tree_type), POINTER                 :: elem
    TYPE(message_send), POINTER              :: m_send
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_scf_step_ener', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(m_send),cp_failure_level,routineP,error,failure)

    elem%scf_energies(MOD(elem%scf_energies_count,4)+1) = m_send%task_real(1)
    elem%scf_energies_count = elem%scf_energies_count +1

  END SUBROUTINE read_scf_step_ener

  !============================================================================
  ! message for broadcasting atom types
  !============================================================================
! *****************************************************************************
!> \brief routines send atom names to the global master 
!>        (using broadcast in a specialized group consisting of the master 
!>        and the first energy worker master)
!> \param atoms
!> \param para_env the communicator environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE communicate_atom_types(atoms, source, para_env, error)
    TYPE(tmc_atom_type), DIMENSION(:), &
      POINTER                                :: atoms
    INTEGER                                  :: source
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'communicate_atom_types', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      ALLOCATABLE, DIMENSION(:)              :: msg(:)
    INTEGER                                  :: i
    LOGICAL                                  :: failure

    CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(source.GE.0,cp_failure_level,routineP,error,failure)
    CPPrecondition(source.LT.para_env%num_pe,cp_failure_level,routineP,error,failure)

    ALLOCATE(msg(SIZE(atoms)))
    IF(para_env%mepos.EQ.source) THEN
      DO i=1, SIZE(atoms)
        msg(i) = atoms(i)%name
      END DO
      CALL mp_bcast(msg, source, para_env%group)
    ELSE
      CALL mp_bcast(msg, source, para_env%group)
      DO i=1, SIZE(atoms)
        atoms(i)%name = msg(i)
      END DO
    END IF
    DEALLOCATE(msg)
  END SUBROUTINE communicate_atom_types

  !============================================================================
  ! message for cancelation or stopping
  !============================================================================
! *****************************************************************************
!> \brief routine cancel the other group participants
!> \param para_env the communicator environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE cancel_group_participants(para_env, tmc_params, error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cancel_group_participants', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: cancel_stat, comm_size, &
                                                participant
    LOGICAL                                  :: failure, flag

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(para_env%group.EQ.0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

    comm_size = para_env%num_pe-1
    !-- send canceling message to all participants
    paritcipant_loop: DO participant=1, comm_size
       cancel_stat = TMC_CANCELING_MESSAGE
       CALL tmc_message(msg_type=cancel_stat, send_recv=send_msg, &
                        dest=participant, para_env=para_env, &
                        tmc_params=tmc_params, error=error)
    END DO paritcipant_loop

    !-- receive canceling receipt from all
    receipt_loop: DO participant=1, comm_size
       CALL tmc_message(msg_type=cancel_stat, send_recv=recv_msg, &
                        dest=participant, para_env=para_env, &
                        tmc_params=tmc_params, &
                        success=flag, wait_for_message=.TRUE., error=error)
       CPPostconditionNoFail(cancel_stat.NE.TMC_CANCELING_RECEIPT,cp_failure_level,routineP,error)
    END DO receipt_loop
  END SUBROUTINE cancel_group_participants

! *****************************************************************************
!> \brief send stop command to all group participants
!> \param 
!> \param 
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 01.2013
! *****************************************************************************
  SUBROUTINE stop_whole_group(para_env, worker_info, tmc_params, error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(elem_array_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: worker_info
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'stop_whole_group', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: dest_rank, stat
    LOGICAL                                  :: failure, flag
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: rank_stoped

!    INTEGER, DIMENSION(MPI_STATUS_SIZE)      :: status_single

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

    ALLOCATE(rank_stoped(0:para_env%num_pe-1))
    rank_stoped(:) = .FALSE.
    rank_stoped(para_env%mepos) = .TRUE.

    ! global master
    IF(PRESENT(worker_info)) THEN
      CPPrecondition(ASSOCIATED(worker_info),cp_failure_level,routineP,error,failure)
      ! canceling running jobs and stop workers
      worker_group_loop: DO dest_rank=1, para_env%num_pe-1
        ! busy workers have to be canceled
        IF(worker_info(dest_rank)%busy) THEN
          stat = TMC_CANCELING_MESSAGE
          CALL tmc_message(msg_type=stat, send_recv=send_msg, dest=dest_rank, &
                           para_env=para_env, tmc_params=tmc_params, &
                           error=error)
        ELSE
        ! send stop message
          stat = TMC_STATUS_FAILED
          CALL tmc_message(msg_type=stat, send_recv=send_msg, dest=dest_rank, &
                           para_env=para_env, tmc_params=tmc_params, &
                           error=error)
        END IF
      END DO worker_group_loop
    ELSE
      ! group master send stop message to all participants
      stat = TMC_STATUS_FAILED
      CALL tmc_message(msg_type=stat, send_recv=send_msg, dest=bcast_group, &
                       para_env=para_env, tmc_params=tmc_params, &
                       error=error)
    END IF

    ! receive stop message receipt
    IF(para_env%mepos.EQ.MASTER_COMM_ID)THEN
      wait_for_receipts: DO
        ! check incomming messages
        stat = TMC_STATUS_WAIT_FOR_NEW_TASK
        dest_rank = 999
        IF(PRESENT(worker_info)) THEN
          ! mast have to be able to receive results, if canceling was too late
          CALL tmc_message(msg_type=stat, send_recv=recv_msg, dest=dest_rank,& 
                           para_env=para_env, tmc_params=tmc_params, &
                           elem_array=worker_info(:), success=flag, error=error)
        ELSE
          CALL tmc_message(msg_type=stat, send_recv=recv_msg, dest=dest_rank, &
                           para_env=para_env, tmc_params=tmc_params, &
                           error=error)
        END IF
        SELECT CASE (stat)
        CASE(TMC_STATUS_WAIT_FOR_NEW_TASK)
          ! no message received
        CASE(TMC_CANCELING_RECEIPT)
          IF(PRESENT(worker_info)) THEN
            worker_info(dest_rank)%busy = .FALSE.
            stat = TMC_STATUS_FAILED
            ! send stop message
            CALL tmc_message(msg_type=stat, send_recv=send_msg, dest=dest_rank,&
                             para_env=para_env, tmc_params=tmc_params, &
                             error=error)
          ELSE
            STOP "group master should not receive cancel receipt"
          END IF
        CASE(TMC_STATUS_STOP_RECEIPT)
          rank_stoped(dest_rank) = .TRUE.
        CASE(TMC_STAT_ENERGY_RESULT, TMC_STAT_NMC_RESULT, TMC_STAT_MD_RESULT, &
             TMC_STAT_SCF_STEP_ENER_RECEIVE, TMC_STAT_APPROX_ENERGY_RESULT, TMC_STAT_ANALYSIS_RESULT)
          ! nothing to do, canceling message already sent
        CASE DEFAULT
          CALL cp_assert(.FALSE., cp_failure_level,cp_assertion_failed,routineP,&
                         "master received status "//cp_to_string(stat)//&
                         " while stopping workers",&
                         error, failure=flag)
        END SELECT
        IF(ALL(rank_stoped)) EXIT wait_for_receipts
      END DO wait_for_receipts
    ELSE
      CALL cp_assert(.FALSE., cp_failure_level,cp_assertion_failed,routineP,&
                     "only (group) master should stop other participants",&
                     error, failure=flag)
    END IF
  END SUBROUTINE stop_whole_group

END MODULE tmc_messages
