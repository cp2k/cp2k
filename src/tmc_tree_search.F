!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief tree nodes search etc. 
!> \par History
!>      11.2012 created [Mandes Schönherr]
!> \author Mandes
! *****************************************************************************

MODULE tmc_tree_search
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE tmc_stati,                       ONLY: TMC_STATUS_WAIT_FOR_NEW_TASK
  USE tmc_tree_references,             ONLY: add_to_references,&
                                             search_and_remove_reference_in_list
  USE tmc_tree_types,                  ONLY: &
       elem_array_type, global_tree_type, status_accepted, &
       status_accepted_result, status_calc_approx_ener, status_calculate_MD, &
       status_calculate_NMC_steps, status_calculate_energy, &
       status_calculated, status_cancel_ener, status_cancel_nmc, &
       status_canceled_ener, status_canceled_nmc, status_created, &
       status_deleted, status_deleted_result, status_rejected, &
       status_rejected_result, tree_type
  USE tmc_types,                       ONLY: tmc_env_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tmc_tree_search'

  PUBLIC :: most_prob_end
  PUBLIC :: search_next_energy_calc
  PUBLIC :: search_canceling_elements
  PUBLIC :: search_parent_element, get_subtree_elements_to_check
  PUBLIC :: search_next_gt_element_to_check
  PUBLIC :: search_end_of_clean_g_tree, search_end_of_clean_tree
  PUBLIC :: count_prepared_nodes_in_trees, count_nodes_in_trees
CONTAINS

  !============================================================================
  ! search tree node
  !============================================================================
! *****************************************************************************
!> \brief search most probable end in global tree to create a new tree node
!>         using the acceptance probabilities for each move type 
!>          of each temperature
!>        routine distinguishes the search for most probable node 
!>         for energy and most probable node with open end 
!>         for new configuration
!>        In case of searching open end: 
!>         routine stops in branch with canceled NMC, 
!>         using this a one possibility
!> \param global_tree_elem starting point for search
!> \param prob return value, the probability of reaching the tree node
!> \param n_acc drection of branch the next tree node should extend
!> \parma search_energy_node flag if configuration for calculating exact 
!>        energy sould be searched
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  RECURSIVE SUBROUTINE most_prob_end(global_tree_elem, prob, n_acc, &
                                     search_energy_node, error)
    TYPE(global_tree_type), POINTER          :: global_tree_elem
    REAL(KIND=dp), INTENT(OUT)               :: prob
    LOGICAL, INTENT(INOUT)                   :: n_acc
    LOGICAL, OPTIONAL                        :: search_energy_node
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'most_prob_end', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: check_accepted, &
                                                check_rejected, failure, &
                                                keep_on, tmp_acc, tmp_nacc
    REAL(KIND=dp)                            :: prob_n_acc, prob_n_nacc
    TYPE(global_tree_type), POINTER          :: ptr_acc, ptr_nacc
    TYPE(tree_type), POINTER                 :: st_elem

    failure = .FALSE.
    NULLIFY(st_elem, ptr_acc, ptr_nacc)

    prob_n_acc  = -100000
    prob_n_nacc = -100000
    check_accepted = .FALSE.
    check_rejected = .FALSE.
    keep_on = .TRUE.

    CPPrecondition(ASSOCIATED(global_tree_elem),cp_failure_level,routineP,error,failure)
    st_elem => global_tree_elem%conf(global_tree_elem%mv_conf)%elem
    CPPrecondition(ASSOCIATED(st_elem),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    !-- follow trajectory until end 
    !-- evaluate following elements using status, and probabilites
    SELECT CASE(global_tree_elem%stat)
    CASE(status_accepted, status_accepted_result)
      check_accepted=.TRUE.
    CASE(status_rejected, status_rejected_result)
      check_rejected=.TRUE.
    CASE DEFAULT
      !-- set directions of searching
      SELECT CASE(st_elem%stat)
      CASE(status_created, status_canceled_ener)
        ! just for searching next element to calculate energy for (found)
        IF(PRESENT(search_energy_node)) THEN
          prob = 0.0_dp ! = log(1) 
          n_acc = .FALSE. ! not needed, but maybe for initialisation
          keep_on = .FALSE.
        ELSE
          check_accepted=.TRUE.
          check_rejected=.TRUE.
        END IF
      CASE(status_canceled_nmc)
        ! just for search new element to create (found)
        ! canceled elements can be reactivated
        ! the parent element is returned, 
        !  the create_new_pt_tree_node check for existing of this node
        IF(.NOT.PRESENT(search_energy_node)) THEN
          prob = 0.0_dp
          n_acc = ASSOCIATED(global_tree_elem%parent%acc, global_tree_elem)
          global_tree_elem => global_tree_elem%parent
          keep_on = .FALSE.
        END IF
      CASE(status_calculated, status_calculate_energy, &
           status_accepted_result, status_accepted, &
           status_rejected, status_rejected_result)
           ! status accepted and rejection needed for swapped 
           !  configurations in parallel tempering
        check_accepted=.TRUE.
        check_rejected=.TRUE.
      CASE(status_calculate_MD, status_calculate_NMC_steps, &
           status_calc_approx_ener)
        ! just for searching next element to create
        IF(.NOT.PRESENT(search_energy_node)) THEN
          check_rejected=.TRUE.
        END IF
      CASE(status_cancel_nmc, status_cancel_ener)
      CASE DEFAULT
        CALL cp_assert(.FALSE.,&
                    cp_failure_level,cp_assertion_failed,routineP,&
                    "unknown sub tree element status "//&
                    cp_to_string(st_elem%stat),error)
      END SELECT
    END SELECT

    IF(keep_on) THEN
      !-- recursive search, remembering lowest element (tree end), 
      !     and multiply probabilities to go there
      !-- search in ACCEPTED branch
      IF(check_accepted)THEN
        ! test if probable accepted child exist and is not rejected
        IF(ASSOCIATED(global_tree_elem%acc))THEN
          ptr_acc => global_tree_elem%acc
          IF(PRESENT(search_energy_node)) THEN
            CALL most_prob_end(global_tree_elem=ptr_acc, prob=prob_n_acc, &
                               n_acc=tmp_acc, &
                               search_energy_node=search_energy_node, &
                               error=error)
          ELSE
            CALL most_prob_end(global_tree_elem=ptr_acc, prob=prob_n_acc, &
                               n_acc=tmp_acc, error=error)
          END IF
          !-- do probability multiplication 
          !    (in logscale because of realy small probabilities)
          prob_n_acc = prob_n_acc + LOG(global_tree_elem%prob_acc)
        ELSE
          ! prob of going in acc or rej direction is 
          !   calculated in parent element
          prob_n_acc = LOG(global_tree_elem%prob_acc)
          IF(PRESENT(search_energy_node)) prob_n_acc = -100000
          ptr_acc => global_tree_elem
          tmp_acc = .TRUE.
        END IF
      END IF

      !-- search in REJECTED branch
      IF(check_rejected)THEN
        ! test if probabliy rejected child exist
        IF(ASSOCIATED(global_tree_elem%nacc))THEN
          ptr_nacc => global_tree_elem%nacc
          IF(PRESENT(search_energy_node)) THEN
            CALL most_prob_end(global_tree_elem=ptr_nacc, prob=prob_n_nacc, &
                               n_acc=tmp_nacc, &
                               search_energy_node=search_energy_node, &
                               error=error)
          ELSE
            CALL most_prob_end(global_tree_elem=ptr_nacc, prob=prob_n_nacc, &
                               n_acc=tmp_nacc, error=error)
          END IF
          !-- do probability multiplication 
          !     (in logscale because of realy small probabilities)
          prob_n_nacc =  prob_n_nacc + LOG(1-global_tree_elem%prob_acc)
        ELSE
          ! prob of going in acc or rej direction is 
          !   calculated in parent element
          prob_n_nacc = LOG(1-global_tree_elem%prob_acc)
          IF(PRESENT(search_energy_node)) prob_n_nacc = -100000
          ptr_nacc => global_tree_elem
          tmp_nacc = .FALSE.
        END IF
      END IF
      ! test wich direction is more likely 
      !   and set result pointer and probability,
      ! remembering the direction 
      IF(prob_n_acc.GE.prob_n_nacc)THEN
        prob = prob_n_acc
        global_tree_elem => ptr_acc
        n_acc = tmp_acc
      ELSE
        prob = prob_n_nacc
        global_tree_elem => ptr_nacc
        n_acc = tmp_nacc
      END IF
    END IF
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE most_prob_end

! *****************************************************************************
!> \brief gt_head head of the global tree
!> \param new_gt_elem return value the energy should be calculated for
!> \param stat routine status return value
!> \param react_count reactivation counter
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE search_next_energy_calc(gt_head, new_gt_elem, stat, react_count, &
                                     error)
    TYPE(global_tree_type), POINTER          :: gt_head, new_gt_elem
    INTEGER                                  :: stat, react_count
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'search_next_energy_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure, flag
    REAL(KIND=dp)                            :: prob

    prob = 0.0_dp
    flag = .FALSE.
    failure = .FALSE.
    CPPrecondition(ASSOCIATED(gt_head),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    new_gt_elem => gt_head

    CALL most_prob_end(global_tree_elem=new_gt_elem, prob=prob, n_acc=flag, &
                       search_energy_node=.TRUE., error=error)

    stat = status_created
    ! set status for master 
    !   (if TMC_STATUS_WAIT_FOR_NEW_TASK, no calculation neccessary)
    IF(.NOT.ASSOCIATED(new_gt_elem).OR.(EXP(prob).LT.1.0E-10)) THEN
      stat=TMC_STATUS_WAIT_FOR_NEW_TASK
    ELSE
      ! reactivate canceled elements
      IF(new_gt_elem%conf(new_gt_elem%mv_conf)%elem%stat.EQ.&
         status_canceled_ener)THEN
        CALL add_to_references(gt_elem=new_gt_elem, error=error)
        react_count = react_count+1
        new_gt_elem%conf(new_gt_elem%mv_conf)%elem%stat = status_created
      END IF
      ! if elem status is not status_created 
      IF(new_gt_elem%conf(new_gt_elem%mv_conf)%elem%stat.NE.status_created)THEN
        stat=TMC_STATUS_WAIT_FOR_NEW_TASK
      END IF
    END IF
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE search_next_energy_calc

! *****************************************************************************
!> \brief searching the parent element (last accepted configuration before)
!> \param current actual tree element 
!> \param parent parent tree element (last accepted one)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  !-- routine searches last (assumed) accepted element in subtree
  RECURSIVE FUNCTION search_parent_element(current, error) RESULT(parent)
    TYPE(tree_type), POINTER                 :: current
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(tree_type), POINTER                 :: parent

    CHARACTER(LEN=*), PARAMETER :: routineN = 'search_parent_element', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(current),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    IF(ASSOCIATED(current%parent))THEN
       ! the result value if the child (we came from) is in acc direction
       parent => current%parent
       IF(ASSOCIATED(parent%nacc, current)) THEN
          parent => search_parent_element(parent, error)
       END IF
    ELSE
       ! if parent not exist, we are at the head of the tree
       parent => current
    END IF
    ! end the timing
    CALL timestop(handle)
    CPPostcondition(ASSOCIATED(parent),cp_failure_level,routineP,error,failure)
  END FUNCTION search_parent_element

! *****************************************************************************
!> \brief search the next global element in the Markov Chain to check
!> \param ptr start point for search, should be on the known Markov Chain
!> \param found flag if routine was successful
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  RECURSIVE SUBROUTINE search_next_gt_element_to_check(ptr, found, error)
    TYPE(global_tree_type), POINTER          :: ptr
    LOGICAL                                  :: found
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'search_next_gt_element_to_check', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    found = .FALSE.

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(ptr),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    ! -- global tree status is not updated after receiving calculations 
    !    (not intrinsically), hence try to check elements with could be ready
    SELECT CASE(ptr%stat)
    CASE(status_accepted_result)
      IF(ASSOCIATED(ptr%acc))THEN
         ptr => ptr%acc
         CALL search_next_gt_element_to_check(ptr, found, error)
      END IF
    CASE(status_rejected_result)
      IF(ASSOCIATED(ptr%nacc))THEN
         ptr => ptr%nacc
         CALL search_next_gt_element_to_check(ptr, found, error)
      END IF
    CASE(status_calculate_energy, status_created, &
        status_calculate_MD,status_calculated,status_calculate_NMC_steps, &
        status_calc_approx_ener, status_accepted, status_rejected)
      found = .TRUE.
    CASE(status_cancel_nmc, status_cancel_ener, status_canceled_nmc, &
         status_canceled_ener)
      ! nothing to do
    CASE DEFAULT
      CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                  routineP,"unexpected status "//cp_to_string(ptr%stat)//&
                  "of global tree elem "//cp_to_string(ptr%nr), &
                  error,failure)
    END SELECT
    ! end the timing
    CALL timestop(handle)

    CPPostcondition(ASSOCIATED(ptr),cp_failure_level,routineP,error,failure)
  END SUBROUTINE search_next_gt_element_to_check

! *****************************************************************************
!> \brief get the changed element of the actual global tree element and its 
!>        related last accepted parent
!> \param gt_act_elem actual global tree element
!> \param elem1, elem2 two subtree elements which should be compared
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE get_subtree_elements_to_check(gt_act_elem, elem1, elem2, error)
    TYPE(global_tree_type), POINTER          :: gt_act_elem
    TYPE(tree_type), INTENT(OUT), POINTER    :: elem1, elem2
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'get_subtree_elements_to_check', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(gt_act_elem),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    IF(gt_act_elem%swaped)THEN
       !------------------------------------------------------------
       !-- take the last accepted configurations for check of both configurations, because 
       !-- in case of swapping, the last accepted elements have to be compared
       IF(gt_act_elem%conf_n_acc(gt_act_elem%conf(gt_act_elem%mv_conf)%elem%sub_tree_nr))THEN
          elem1 => gt_act_elem%conf(gt_act_elem%mv_conf)%elem
       ELSE
          elem1 => search_parent_element(gt_act_elem%conf(gt_act_elem%mv_conf)%elem, error)
       END IF
       ! second element
       IF(gt_act_elem%conf_n_acc(gt_act_elem%conf(gt_act_elem%mv_conf+1)%elem%sub_tree_nr))THEN
          elem2 => gt_act_elem%conf(gt_act_elem%mv_conf+1)%elem
       ELSE
          elem2 => search_parent_element(gt_act_elem%conf(gt_act_elem%mv_conf+1)%elem, error)
       END IF
    ELSE
       elem1 => gt_act_elem%conf(gt_act_elem%mv_conf)%elem
       elem2 => search_parent_element(elem1, error)
    END IF

    ! end the timing
    CALL timestop(handle)

    CPPostcondition(ASSOCIATED(gt_act_elem),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(elem1),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(elem2),cp_failure_level,routineP,error,failure)
  END SUBROUTINE get_subtree_elements_to_check

! *****************************************************************************
!> \brief searches last element on trajectory,
!>        until where the sides of the tree are deleted (of global tree)
!>        also found the last accepted element before
!> \param last_acc returns last accepted element in cleaned tree part
!> \param tree_ptr end point of search
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  RECURSIVE SUBROUTINE search_end_of_clean_g_tree(last_acc, tree_ptr, error)
    TYPE(global_tree_type), POINTER          :: last_acc, tree_ptr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'search_end_of_clean_g_tree', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(last_acc),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tree_ptr),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    SELECT CASE(tree_ptr%stat)
    CASE(status_accepted_result)
      IF(ASSOCIATED(tree_ptr%acc).AND..NOT.ASSOCIATED(tree_ptr%nacc))THEN
        last_acc => tree_ptr
        tree_ptr => tree_ptr%acc
        CALL search_end_of_clean_g_tree(last_acc, tree_ptr, error)
      END IF
    CASE(status_rejected_result)
      IF(ASSOCIATED(tree_ptr%nacc).AND..NOT.ASSOCIATED(tree_ptr%acc))THEN
        tree_ptr => tree_ptr%nacc
        CALL search_end_of_clean_g_tree(last_acc, tree_ptr, error)
      END IF
    CASE(status_calculated, status_calculate_energy, status_created, status_accepted, status_rejected, &
         status_calculate_MD, status_calculate_NMC_steps, status_calc_approx_ener, &
         status_canceled_ener, status_canceled_nmc, &
         status_cancel_nmc, status_cancel_ener)
       ! nothing to do
    CASE DEFAULT
      CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                routineP,"the global tree element "//cp_to_string(tree_ptr%nr)//&
                " stat "//cp_to_string(tree_ptr%stat)//" is UNknown",&
                error,failure)
    END SELECT
    ! end the timing
    CALL timestop(handle)
    CPPostcondition(ASSOCIATED(last_acc),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(tree_ptr),cp_failure_level,routineP,error,failure)
  END SUBROUTINE search_end_of_clean_g_tree

! *****************************************************************************
!> \brief searches last element on trajectory,
!>        until where the sides of the tree are deleted (in sub tree)
!>        also found the last accepted element before.
!>        searches the last element which at least have ONE (not calculated) 
!>        node in the tree branch
!> \param tree_ptr 
!> \param last_acc
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  RECURSIVE SUBROUTINE search_end_of_clean_tree(tree_ptr, last_acc, error)
    TYPE(tree_type), POINTER                 :: tree_ptr
    TYPE(tree_type), INTENT(IN), POINTER     :: last_acc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'search_end_of_clean_tree', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(tree_ptr),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(last_acc),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    IF(.NOT.ASSOCIATED(last_acc, tree_ptr)) THEN
       IF(ASSOCIATED(tree_ptr%acc).AND..NOT.ASSOCIATED(tree_ptr%nacc))THEN
         tree_ptr => tree_ptr%acc
         CALL search_end_of_clean_tree(tree_ptr, last_acc, error)
       ELSE IF(ASSOCIATED(tree_ptr%nacc).AND..NOT.ASSOCIATED(tree_ptr%acc))THEN
         tree_ptr => tree_ptr%nacc
         CALL search_end_of_clean_tree(tree_ptr, last_acc, error)
       END IF
    END IF
    ! end the timing
    CALL timestop(handle)
    CPPostcondition(ASSOCIATED(tree_ptr),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(last_acc),cp_failure_level,routineP,error,failure)
  END SUBROUTINE search_end_of_clean_tree

! *****************************************************************************
!> \brief searches in all branches down below the entered global tree element 
!>        for elements to cancel, if prob is present start searching at the 
!>        related tree child node
!> \param pt_elem_in start search point 
!> \param prob the acceptance probability of the tree element to define 
!>        the direction to start with
!> \param tmc_env TMC environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  RECURSIVE SUBROUTINE search_canceling_elements(pt_elem_in, prob, tmc_env,&
                                                 error)
    TYPE(global_tree_type), INTENT(IN), &
      POINTER                                :: pt_elem_in
    REAL(KIND=dp), OPTIONAL                  :: prob
    TYPE(tmc_env_type), POINTER              :: tmc_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'search_canceling_elements', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure, ready
    TYPE(global_tree_type), POINTER          :: act_pt_ptr, pt_elem

    failure = .FALSE.
    NULLIFY(pt_elem, act_pt_ptr)
    CPPrecondition(ASSOCIATED(pt_elem_in),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_env),cp_failure_level,routineP,error,failure)
    
    ! start the timing
    CALL timeset(routineN,handle)

    ready = .TRUE.
    ! if prob present select the related branch
    IF(PRESENT(prob)) THEN
       IF(prob.LT.1.0E-10 .AND.ASSOCIATED(pt_elem_in%acc)) THEN
          pt_elem => pt_elem_in%acc
       ELSE IF(prob.GT.(1.0_dp-1.0E-10) .AND.ASSOCIATED(pt_elem_in%nacc)) THEN
          pt_elem => pt_elem_in%nacc
       ELSE
         ready = .FALSE.
       END IF
    ELSE
       pt_elem => pt_elem_in
    END IF

    IF(ready) THEN
      IF(ASSOCIATED(pt_elem%conf(pt_elem%mv_conf)%elem))THEN
         SELECT CASE(pt_elem%conf(pt_elem%mv_conf)%elem%stat)
         CASE(status_accepted_result, status_accepted, status_rejected_result,&
              status_rejected, status_created, status_cancel_nmc, &
              status_cancel_ener, status_canceled_nmc, status_canceled_ener, &
              status_calculated, status_deleted, status_deleted_result, &
              status_calc_approx_ener) ! no canceling
         CASE(status_calculate_NMC_steps, status_calculate_MD, &
              status_calculate_energy)
           CALL search_and_remove_reference_in_list(gt_ptr=pt_elem, &
                   elem=pt_elem%conf(pt_elem%mv_conf)%elem, tmc_env=tmc_env, &
                   error=error)

         CASE DEFAULT
           CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                   routineP,"unknown status of subtree element"//&
                   cp_to_string(pt_elem%conf(pt_elem%mv_conf)%elem%stat),&
                   error,failure)
         END SELECT
      END IF
      !-- go until the ends ot he tree, to search for elements to cancel
      !-- check if child nodes exist
      IF(ASSOCIATED(pt_elem%acc))THEN
         act_pt_ptr => pt_elem%acc
         CALL search_canceling_elements(pt_elem_in=act_pt_ptr,tmc_env=tmc_env,&
                                        error=error)
      END IF
      IF(ASSOCIATED(pt_elem%nacc))THEN
         act_pt_ptr => pt_elem%nacc
         CALL search_canceling_elements(pt_elem_in=act_pt_ptr,tmc_env=tmc_env,&
                                        error=error)
      END IF
    END IF
    ! end the timing
    CALL timestop(handle)
    CPPostcondition(ASSOCIATED(pt_elem_in),cp_failure_level,routineP,error,failure)
  END SUBROUTINE search_canceling_elements

! *****************************************************************************
!> \brief searches for created configurations in all subtrees
!> \param global_tree_ptr pointer to one global tree element
!> \param counters array returning the counters for each subtree
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 01.2013
! *****************************************************************************
  SUBROUTINE count_prepared_nodes_in_trees(global_tree_ptr, counters, error)
    TYPE(global_tree_type), INTENT(IN), &
      POINTER                                :: global_tree_ptr
    INTEGER, DIMENSION(:), INTENT(OUT), &
      POINTER                                :: counters
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'count_prepared_nodes_in_trees', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    LOGICAL                                  :: failure
    TYPE(tree_type), POINTER                 :: begin_ptr

    failure = .FALSE.
    NULLIFY(begin_ptr)

    CPPrecondition(ASSOCIATED(global_tree_ptr),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(counters),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(counters(1:)).EQ.SIZE(global_tree_ptr%conf(:)),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    counters(:) = 0
    DO i=1, SIZE(global_tree_ptr%conf(:))
       begin_ptr => global_tree_ptr%conf(i)%elem
       CALL count_prepared_nodes_in_subtree(tree_ptr=begin_ptr, &
                                            counter=counters(i), error=error)
    END DO

    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE count_prepared_nodes_in_trees

! *****************************************************************************
!> \brief counts the prepared tree nodes in subtrees
!> \param tree_ptr pointer to one subtree element 
!> \param counter returning the amount of prepared 
!>        (ready for energy calculation) elements ind certain sub tree
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 01.2013
! *****************************************************************************
  RECURSIVE SUBROUTINE count_prepared_nodes_in_subtree(tree_ptr, counter,error)
    TYPE(tree_type), POINTER                 :: tree_ptr
    INTEGER                                  :: counter
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'count_prepared_nodes_in_subtree', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(tree_type), POINTER                 :: tmp_ptr

    failure = .FALSE.
    NULLIFY(tmp_ptr)

    CPPrecondition(ASSOCIATED(tree_ptr),cp_failure_level,routineP,error,failure)

    SELECT CASE(tree_ptr%stat)
    CASE(status_accepted, status_accepted_result)
      IF(ASSOCIATED(tree_ptr%acc)) THEN
        tmp_ptr => tree_ptr%acc
        CALL count_prepared_nodes_in_subtree(tmp_ptr, counter, error)
      END IF
    CASE(status_rejected, status_rejected_result)
      IF(ASSOCIATED(tree_ptr%nacc)) THEN
        tmp_ptr => tree_ptr%nacc
        CALL count_prepared_nodes_in_subtree(tmp_ptr, counter, error)
      END IF
    CASE(status_created, status_calculate_MD, status_calculate_NMC_steps, &
         status_calc_approx_ener, status_calculated, status_calculate_energy)
      IF(tree_ptr%stat.EQ.status_created) counter = counter + 1
      IF(ASSOCIATED(tree_ptr%acc)) THEN
        tmp_ptr => tree_ptr%acc
        CALL count_prepared_nodes_in_subtree(tmp_ptr, counter, error)
      END IF
      IF(ASSOCIATED(tree_ptr%nacc)) THEN
        tmp_ptr => tree_ptr%nacc
        CALL count_prepared_nodes_in_subtree(tmp_ptr, counter, error)
      END IF
    CASE(status_cancel_nmc, status_cancel_ener, status_canceled_nmc, &
         status_canceled_ener)
      !TODO maybe also count caneled confs for debug output
    CASE DEFAULT
      CALL cp_assert(.FALSE., cp_failure_level,cp_assertion_failed,routineP,&
                   "stat "//cp_to_string(tree_ptr%stat)//&
                   "of elem "//cp_to_string(tree_ptr%nr)//&
                   "unknown.",&
                   error)
    END SELECT
  END SUBROUTINE count_prepared_nodes_in_subtree

! *****************************************************************************
!> \brief counts the number of existing nodes in global and subtrees
!> \param global_tree_ptr pointer to one global tree element
!> \param end_of_clean_trees points to the last elements of the clean sub trees
!> \param counters array returning the counters for each subtree
!> \param head_elements_nr node number of the existing 
!>        global and sub tree heads
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 01.2013
! *****************************************************************************
  SUBROUTINE count_nodes_in_trees(global_tree_ptr, end_of_clean_trees, &
                                  counters, head_elements_nr, error)
    TYPE(global_tree_type), POINTER          :: global_tree_ptr
    TYPE(elem_array_type), DIMENSION(:), &
      POINTER                                :: end_of_clean_trees
    INTEGER, DIMENSION(:), INTENT(OUT), &
      POINTER                                :: counters, head_elements_nr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'count_nodes_in_trees', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    LOGICAL                                  :: failure
    TYPE(global_tree_type), POINTER          :: begin_gt_ptr
    TYPE(tree_type), POINTER                 :: begin_ptr

    failure = .FALSE.
    NULLIFY(begin_gt_ptr, begin_ptr)

    CPPrecondition(ASSOCIATED(global_tree_ptr),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(end_of_clean_trees),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(counters),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(counters(1:)).EQ.SIZE(global_tree_ptr%conf(:)),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    begin_gt_ptr => global_tree_ptr
    counters(:) = 0
    DO
      IF(.NOT.ASSOCIATED(begin_gt_ptr%parent)) EXIT
      begin_gt_ptr => begin_gt_ptr%parent
    END DO
    head_elements_nr(0) = begin_gt_ptr%nr
    CALL count_nodes_in_global_tree(begin_gt_ptr, counters(0), error)
    DO i=1, SIZE(end_of_clean_trees(:))
      begin_ptr => end_of_clean_trees(i)%elem
      DO
        IF(.NOT.ASSOCIATED(begin_ptr%parent)) EXIT
        begin_ptr => begin_ptr%parent
      END DO
      head_elements_nr(i) = begin_ptr%nr
      CALL count_nodes_in_tree(begin_ptr, counters(i), error)
    END DO

    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE count_nodes_in_trees

! *****************************************************************************
!> \brief counts existing nodes in global tree
!> \param ptr global tree head 
!> \param counter return value with the amount of existing global tree elements
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 01.2013
! *****************************************************************************
  RECURSIVE SUBROUTINE count_nodes_in_global_tree(ptr, counter, error)
    TYPE(global_tree_type), INTENT(IN), &
      POINTER                                :: ptr
    INTEGER, INTENT(INOUT)                   :: counter
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'count_nodes_in_global_tree', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(ptr),cp_failure_level,routineP,error,failure)

    counter = counter +1

    IF(ASSOCIATED(ptr%acc)) &
      CALL count_nodes_in_global_tree(ptr%acc, counter, error)
    IF(ASSOCIATED(ptr%nacc)) &
      CALL count_nodes_in_global_tree(ptr%nacc, counter, error)
  END SUBROUTINE count_nodes_in_global_tree

! *****************************************************************************
!> \brief counts existing nodes in certain sub tree
!> \param ptr subtree tree head 
!> \param counter return value with the amount of existing sub tree elements
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 01.2013
! *****************************************************************************
  RECURSIVE SUBROUTINE count_nodes_in_tree(ptr, counter, error)
    TYPE(tree_type), POINTER                 :: ptr
    INTEGER                                  :: counter
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'count_nodes_in_tree', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(ptr),cp_failure_level,routineP,error,failure)

    counter = counter +1

    IF(ASSOCIATED(ptr%acc)) &
      CALL count_nodes_in_tree(ptr%acc, counter, error)
    IF(ASSOCIATED(ptr%nacc)) &
      CALL count_nodes_in_tree(ptr%nacc, counter, error)
  END SUBROUTINE count_nodes_in_tree
END MODULE tmc_tree_search
