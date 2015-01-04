!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief global tree references
!>        - BECAUSE acceptance check use global tree randon numbers and
!>            (in case of parallel tempering) several global tree node refer to a 
!>            single sub tree node (which is the changed one in the global tree)
!>        - the references are used to update the global tree acceptance probability
!>            for every global tree element seperately
!>        Hence a list of all global tree nodes, using the related subtree node, 
!>            is created.
!> \par History
!>      11.2012 created [Mandes Schönherr]
!> \author Mandes
! *****************************************************************************

MODULE tmc_tree_references
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE tmc_cancelation,                 ONLY: add_to_canceling_list
  USE tmc_tree_types,                  ONLY: global_tree_type,&
                                             gt_elem_list_type,&
                                             tree_type
  USE tmc_types,                       ONLY: tmc_env_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tmc_tree_references'

  PUBLIC :: add_to_references
  PUBLIC :: search_and_remove_reference_in_list
  PUBLIC :: remove_subtree_element_of_all_references
  PUBLIC :: remove_gt_references
CONTAINS

! *****************************************************************************
!> \brief adds global tree reference to the modified sub tree element(s)
!> \param gt_elem actual global tree element
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE add_to_references(gt_elem, error)
    TYPE(global_tree_type), POINTER          :: gt_elem
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'add_to_references', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(gt_elem_list_type), POINTER         :: tmp_pt_list_elem

    failure = .FALSE.
    NULLIFY(tmp_pt_list_elem)

    CPPrecondition(ASSOCIATED(gt_elem),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    ! create reference and add at the beginning of the list
    ALLOCATE(tmp_pt_list_elem)
    tmp_pt_list_elem%gt_elem => gt_elem
    IF(ASSOCIATED(gt_elem%conf(gt_elem%mv_conf)%elem%gt_nodes_references)) THEN
      tmp_pt_list_elem%next => gt_elem%conf(gt_elem%mv_conf)%elem%gt_nodes_references
    ELSE
      tmp_pt_list_elem%next => NULL()
    END IF
    gt_elem%conf(gt_elem%mv_conf)%elem%gt_nodes_references => tmp_pt_list_elem

    ! in case of swapped configurations both are necessary to do acceptance probability update 
    !   also when second configuration returns a value
    IF(gt_elem%swaped) THEN
      ! add reference to swapped elem
      ALLOCATE(tmp_pt_list_elem)
      tmp_pt_list_elem%gt_elem => gt_elem
      IF(ASSOCIATED(gt_elem%conf(gt_elem%mv_conf+1)%elem%gt_nodes_references)) THEN
        tmp_pt_list_elem%next => gt_elem%conf(gt_elem%mv_conf+1)%elem%gt_nodes_references
      ELSE
        tmp_pt_list_elem%next => NULL()
      END IF
      gt_elem%conf(gt_elem%mv_conf+1)%elem%gt_nodes_references => tmp_pt_list_elem
    END IF
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE add_to_references

! *****************************************************************************
!> \brief removes the global tree references of this actual global tree element
!>        from all related sub tree elements
!> \param gt_ptr actual global tree element
!> \param tmc_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE remove_gt_references(gt_ptr, tmc_env, error)
    TYPE(global_tree_type), POINTER          :: gt_ptr
    TYPE(tmc_env_type), POINTER              :: tmc_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'remove_gt_references', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(gt_ptr),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_env),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    CALL search_and_remove_reference_in_list(gt_ptr=gt_ptr, &
            elem=gt_ptr%conf(gt_ptr%mv_conf)%elem, tmc_env=tmc_env,error=error)

    ! in case of parallel tempering also the reference in the second swaped configuration has to be removed
    IF(gt_ptr%swaped) THEN
       CALL search_and_remove_reference_in_list(gt_ptr=gt_ptr, &
               elem=gt_ptr%conf(gt_ptr%mv_conf+1)%elem, tmc_env=tmc_env, error=error)
    END IF
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE remove_gt_references

! *****************************************************************************
!> \brief removes the pointers to a certain subtree element from every related
!>        global tree element
!> \param ptr sub tree element
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE remove_subtree_element_of_all_references(ptr, error)
    TYPE(tree_type), POINTER                 :: ptr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'remove_subtree_element_of_all_references', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=2000)                      :: list_of_nr
    INTEGER                                  :: handle, i
    LOGICAL                                  :: failure
    TYPE(gt_elem_list_type), POINTER         :: tmp_gt_list_ptr

    failure = .FALSE.
    NULLIFY(tmp_gt_list_ptr)

    CPPrecondition(ASSOCIATED(ptr),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    pt_node_ref_loop: DO WHILE(ASSOCIATED(ptr%gt_nodes_references))
      tmp_gt_list_ptr => ptr%gt_nodes_references
      CPPrecondition(ASSOCIATED(tmp_gt_list_ptr%gt_elem),cp_failure_level,routineP,error,failure)
      CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                   routineP,"found reference of global tree node "//&
                   cp_to_string(tmp_gt_list_ptr%gt_elem%nr)//&
                   ", while removing sub tree node "//&
                   cp_to_string(ptr%sub_tree_nr)//cp_to_string(ptr%nr),&
                   error,failure)
      ! check if configurations exist
      IF(ASSOCIATED(tmp_gt_list_ptr%gt_elem%conf(tmp_gt_list_ptr%gt_elem%mv_conf)%elem)) THEN
        IF(ASSOCIATED(ptr, tmp_gt_list_ptr%gt_elem%conf(tmp_gt_list_ptr%gt_elem%mv_conf)%elem)) THEN
          tmp_gt_list_ptr%gt_elem%conf(tmp_gt_list_ptr%gt_elem%mv_conf)%elem => NULL()
        ! in case of swapping the second configuration could be the related one
        ELSE IF(ASSOCIATED(ptr, tmp_gt_list_ptr%gt_elem%conf(tmp_gt_list_ptr%gt_elem%mv_conf+1)%elem)) THEN
          tmp_gt_list_ptr%gt_elem%conf(tmp_gt_list_ptr%gt_elem%mv_conf+1)%elem => NULL()
        ELSE
          list_of_nr = ""
          DO i=1, SIZE(tmp_gt_list_ptr%gt_elem%conf)
            WRITE(list_of_nr,*)TRIM(ADJUSTL(list_of_nr)), tmp_gt_list_ptr%gt_elem%conf(i)%elem%sub_tree_nr, &
                  tmp_gt_list_ptr%gt_elem%conf(i)%elem%nr, " | "
          END DO
          CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,&
                   routineP,"for subtree "//&
                   cp_to_string(ptr%sub_tree_nr)//&
                   "element "//cp_to_string(ptr%nr)//&
                   "global tree element"//cp_to_string(tmp_gt_list_ptr%gt_elem%nr)//&
                   "swaped"//cp_to_string(tmp_gt_list_ptr%gt_elem%swaped)//&
                   "moved elem"//cp_to_string(tmp_gt_list_ptr%gt_elem%mv_conf)//&
                   "with the related subtree, elements: "//&
                   TRIM(ADJUSTL(list_of_nr)), &
                   error,failure)
        END IF
      ELSE
        CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,&
                   routineP,"for subtree "//cp_to_string(ptr%sub_tree_nr)//&
                   "element "//cp_to_string(ptr%nr)//&
                   " is not related to global tree node "//cp_to_string(tmp_gt_list_ptr%gt_elem%nr)//&
                   "(anymore).",error,failure)
      END IF
      ptr%gt_nodes_references => ptr%gt_nodes_references%next
      DEALLOCATE(tmp_gt_list_ptr)
      NULLIFY(tmp_gt_list_ptr)
    END DO pt_node_ref_loop

    ! end the timing
    CALL timestop(handle)

    CPPostcondition(.NOT.ASSOCIATED(ptr%gt_nodes_references),cp_failure_level,routineP,error,failure)
  END SUBROUTINE remove_subtree_element_of_all_references

! *****************************************************************************
!> \brief removes the global tree references of this actual global tree element
!>        from all related sub tree elements
!> \param gt_ptr actual global tree element
!> \param elem ...
!> \param tmc_env TMC environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE search_and_remove_reference_in_list(gt_ptr, elem, tmc_env, error)
    TYPE(global_tree_type), POINTER          :: gt_ptr
    TYPE(tree_type), POINTER                 :: elem
    TYPE(tmc_env_type), POINTER              :: tmc_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'search_and_remove_reference_in_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(gt_elem_list_type), POINTER         :: tmp_gt_list_last_ptr, &
                                                tmp_gt_list_ptr

    failure = .FALSE.
    NULLIFY(tmp_gt_list_ptr, tmp_gt_list_last_ptr)

    ! nothing to do, when subtree element is already deleted
    IF(.NOT.ASSOCIATED(elem)) RETURN
    IF(.NOT.ASSOCIATED(gt_ptr)) RETURN

    CPPrecondition(ASSOCIATED(tmc_env),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    ! set the entry point od the list
    tmp_gt_list_ptr => elem%gt_nodes_references
    tmp_gt_list_last_ptr => elem%gt_nodes_references

    ! search related reference
    DO WHILE(ASSOCIATED(tmp_gt_list_ptr))
       ! remove reference, if it is related to the global tree element
       IF(ASSOCIATED(tmp_gt_list_ptr%gt_elem, gt_ptr)) THEN
          ! first reference?
          IF(ASSOCIATED(tmp_gt_list_ptr, elem%gt_nodes_references)) THEN
             ! additionally last reference (the only one)?
             IF(.NOT.ASSOCIATED(tmp_gt_list_ptr%next)) THEN
                ! last element in list -> cancel calculation
                CALL add_to_canceling_list(elem=elem, tmc_env=tmc_env, error=error)
                elem%gt_nodes_references => NULL()
                tmp_gt_list_last_ptr => NULL()
             ELSE
                ! if first list element and NOT last one:
                ! set list pointer to second element
                elem%gt_nodes_references => tmp_gt_list_ptr%next
                tmp_gt_list_last_ptr => elem%gt_nodes_references
             END IF
          ELSE
             ! if NOT first one
             ! skip that element in list
             tmp_gt_list_last_ptr%next => tmp_gt_list_ptr%next
          END IF

          ! deallocate list element
          DEALLOCATE(tmp_gt_list_ptr)
          NULLIFY(tmp_gt_list_ptr)
          ! going back to last list element
          tmp_gt_list_ptr => tmp_gt_list_last_ptr
       END IF
       ! setting to next list element
       tmp_gt_list_last_ptr => tmp_gt_list_ptr
       ! go to next list element, if defined
       IF(ASSOCIATED(tmp_gt_list_ptr)) tmp_gt_list_ptr => tmp_gt_list_ptr%next
    END DO
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE search_and_remove_reference_in_list

END MODULE tmc_tree_references
