!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief module for printing tree structures in GraphViz dot files
!>        for visualizing the trees
!> \par History
!>      12.2012 created [Mandes Schönherr]
!> \author Mandes
! *****************************************************************************
!----------------------------------------------------------------------!
! Tree Monte Carlo (TMC) a program for parallel Monte Carlo simulation
! \author Mandes Schoenherr
!----------------------------------------------------------------------!
MODULE tmc_dot_tree
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE tmc_file_io,                     ONLY: expand_file_name_char,&
                                             expand_file_name_temp
  USE tmc_move_types,                  ONLY: mv_type_swap_conf
  USE tmc_tree_types,                  ONLY: &
       global_tree_type, gt_elem_list_type, status_accepted, &
       status_accepted_result, status_calc_approx_ener, status_calculate_MD, &
       status_calculate_NMC_steps, status_calculate_energy, &
       status_calculated, status_cancel_ener, status_cancel_nmc, &
       status_canceled_ener, status_canceled_nmc, status_created, &
       status_deleted, status_deleted_result, status_rejected, &
       status_rejected_result, tree_type
  USE tmc_types,                       ONLY: tmc_param_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tmc_dot_tree'

  PUBLIC :: init_draw_trees, finalize_draw_tree
  PUBLIC :: create_dot_color, create_global_tree_dot_color
  PUBLIC :: create_dot, create_global_tree_dot

  INTEGER :: DEBUG = 0
!  CHARACTER(LEN=30) :: filename ="tree.dot"

CONTAINS
! *****************************************************************************
!> \brief returns extended filename for global and sub trees
!> \param tmc_params param environment for creating the file name
!> \param ind index of the subtree (0 = global tree)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval filename ...
!> \author Mandes 12.2012
! *****************************************************************************
  FUNCTION get_dot_file_name(tmc_params, ind, error) RESULT(filename)
    TYPE(tmc_param_type), POINTER            :: tmc_params
    INTEGER                                  :: ind
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(LEN=50)                        :: filename

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_dot_file_name', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    filename = ""

    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)
    CPPrecondition(ind.GE.0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params%Temp),cp_failure_level,routineP,error,failure)
    CPPrecondition(ind.LE.SIZE(tmc_params%Temp),cp_failure_level,routineP,error,failure)

    IF(ind.EQ.0) THEN
      filename = TRIM(expand_file_name_char(tmc_params%dot_file_name,"global",error))
    ELSE
      filename = TRIM(expand_file_name_temp(file_name=tmc_params%dot_file_name,&
                                            rvalue=tmc_params%Temp(ind),error=error))
    END IF

    CPPostcondition(filename.NE."",cp_failure_level,routineP,error,failure)
  END FUNCTION get_dot_file_name
! *****************************************************************************
!> \brief initializes the dot files (open and write headers)
!> \param tmc_params param environment for creating the file name
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE init_draw_trees(tmc_params, error)
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_draw_trees', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: file_ptr, i
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

    ! global tree
    CALL open_file(file_name=get_dot_file_name(tmc_params, 0, error), file_status="REPLACE", &
                   file_action="WRITE", unit_number=file_ptr)
    WRITE(file_ptr,*)"digraph G {"
    WRITE(file_ptr,*)'  size="8.27,11.69"'
    CALL write_legend(file_ptr, error)
    CALL close_file(unit_number=file_ptr, keep_preconnection=.TRUE.)
 
    ! subtrees
    DO i=1, SIZE(tmc_params%Temp)
       CALL open_file(file_name=get_dot_file_name(tmc_params, i, error), file_status="REPLACE", &
                   file_action="WRITE", unit_number=file_ptr)
       WRITE(file_ptr,*)"digraph G {"
       WRITE(file_ptr,*)'  size="8.27,11.69"'
       CALL write_legend(file_ptr, error)
       CALL close_file(unit_number=file_ptr, keep_preconnection=.TRUE.)
    END DO
  END SUBROUTINE init_draw_trees

! *****************************************************************************
!> \brief close the dot files (write tails)
!> \param tmc_params param environment for creating the file name
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE finalize_draw_tree(tmc_params, error)
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'finalize_draw_tree', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: file_ptr, i
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

    ! global tree
    CALL open_file(file_name=get_dot_file_name(tmc_params, 0, error), &
                   file_status="OLD", file_action="WRITE", &
                   file_position="APPEND", unit_number=file_ptr)
    WRITE(file_ptr,*)"}"
    CALL close_file(unit_number=file_ptr)

    DO i=1, SIZE(tmc_params%Temp)
      CALL open_file(file_name=get_dot_file_name(tmc_params, i, error), file_status="OLD", &
                     file_action="WRITE", file_position="APPEND", unit_number=file_ptr)
      WRITE(file_ptr,*)"}"
      CALL close_file(unit_number=file_ptr)
    END DO
  END SUBROUTINE finalize_draw_tree

! *****************************************************************************
!> \brief writes the legend in the file
!> \param file_ptr file pointer 
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE write_legend(file_ptr, error)
    INTEGER, INTENT(IN)                      :: file_ptr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_legend', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(file_ptr.GT.0,cp_failure_level,routineP,error,failure)

    WRITE(file_ptr,*)'//LEGEND'
    WRITE(file_ptr,*)'subgraph clusterLegend {'
    WRITE(file_ptr,*)'  label="Legend:" labelloc=t fontsize=30'
    WRITE(file_ptr,*)'  centered=false'
    WRITE(file_ptr,*)'  color=black'
    WRITE(file_ptr,*)'  leg1 -> leg2 -> leg2_2 -> leg2_3 -> leg2_4 -> leg3 -> '//&
                        'leg4 -> leg5 -> leg6 -> leg7_1 -> leg7 -> '// &
                        'leg8_1 -> leg8 -> leg9 -> leg10 [style=invis]'
    WRITE(file_ptr,*)'  {rank=same leg1 [fontsize=30, label="node created"          , color=black]}'
    WRITE(file_ptr,*)'  {rank=same leg2 [fontsize=30, label="configuration created" , style=filled,    color=gray]}'
    WRITE(file_ptr,*)'  {rank=same leg2_2 [fontsize=30, label="calc energy" , style=filled,    color=brown]}'
    WRITE(file_ptr,*)'  {rank=same leg2_2 [fontsize=30, label="calc energy" , style=filled,    color=wheat]}'
    WRITE(file_ptr,*)'  {rank=same leg2_3 [fontsize=30, label="calc HMC" , style=filled,    color=goldenrod]}'
    WRITE(file_ptr,*)'  {rank=same leg2_4 [fontsize=30, label="calc NMC" , style=filled,    color=peru]}'
    WRITE(file_ptr,*)'  {rank=same leg3 [fontsize=30, label="accepted"              , color=greenyellow]}'
    WRITE(file_ptr,*)'  {rank=same leg4 [fontsize=30, label="rejected"              , color=red]}'
    WRITE(file_ptr,*)'  {rank=same leg5 [fontsize=30, label="trajec"                , '//&
                         'style=filled,    color=gold, shape=polygon, sides=4]}'
    WRITE(file_ptr,*)'  {rank=same leg6 [fontsize=30, label="energy calculated"     , '// &
                         'style=filled,    color=blue, fontcolor=white]}'
    WRITE(file_ptr,*)'  {rank=same leg7_1 [fontsize=30, label="cancel NMC send"     , '// &
                         'style=filled,    color=deeppink, fontcolor=white]}'
    WRITE(file_ptr,*)'  {rank=same leg7 [fontsize=30, label="canceled NMC"          , '//&
                        'style=filled,    color=darkorchid1, fontcolor=white]}'
    WRITE(file_ptr,*)'  {rank=same leg8_1 [fontsize=30, label="cancel ENERGY send"    , '// &
                        'style=filled,    color=cornflowerblue]}'
    WRITE(file_ptr,*)'  {rank=same leg8 [fontsize=30, label="canceled ENERGY"       , '// &
                        'style=filled,    color=cyan]}'
    WRITE(file_ptr,*)'  {rank=same leg9 [fontsize=30, label="deleted"               , '// &
                        'style=filled,    shape=polygon, sides=3, color=black,fontcolor=white]}'
    WRITE(file_ptr,*)'  {rank=same leg10 [fontsize=30, label="deleted trajectory"   , '// &
                        'style=filled,    shape=polygon, sides=5, color=gold]}'
    WRITE(file_ptr,*)' }'
  END SUBROUTINE write_legend

! *****************************************************************************
!> \brief write/change color related to certain tree element status 
!> \param node_nr the index of the tree node
!> \param stat tree element status
!> \param filename the filename for the grapgviz dot files
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE write_color(node_nr, stat, filename, error)
    INTEGER                                  :: node_nr, stat
    CHARACTER(LEN=50)                        :: filename
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_color', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=11)                        :: label
    INTEGER                                  :: file_ptr
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(filename.NE."",cp_failure_level,routineP,error,failure)
    CPPrecondition(node_nr.GE.0,cp_failure_level,routineP,error,failure)

    CALL open_file(file_name=filename, file_status="OLD", &
                   file_action="WRITE", file_position="APPEND", unit_number=file_ptr)
    WRITE(label, FMT='(I10,A)')node_nr, "["
    SELECT CASE(stat)
    CASE(status_created)
      WRITE(file_ptr,*) TRIM(label),'style=filled, color=gray]'
    CASE(status_accepted)
      WRITE(file_ptr,*) TRIM(label),'color=green]'
    CASE(status_rejected)
      WRITE(file_ptr,*) TRIM(label),'color=red]'
    CASE(status_accepted_result)
      WRITE(file_ptr,*) TRIM(label),'style=filled, color=green, shape=polygon, sides=4]'
    CASE(status_rejected_result)
      WRITE(file_ptr,*) TRIM(label),'style=filled, color=red, shape=polygon, sides=4]'
    CASE(status_calculated)
      WRITE(file_ptr,*) TRIM(label),'style=filled, color=blue]'
    CASE(status_cancel_nmc)
      WRITE(file_ptr,*) TRIM(label),'style=filled, color=deeppink]'
    CASE(status_cancel_ener)
      WRITE(file_ptr,*) TRIM(label),'style=filled, color=cornflowerblue]'
    CASE(status_canceled_nmc)
      WRITE(file_ptr,*) TRIM(label),'style=filled, color=darkorchid1]'
    CASE(status_canceled_ener)
      WRITE(file_ptr,*) TRIM(label),'style=filled, color=cyan]'
    CASE(status_deleted)
      WRITE(file_ptr,*) TRIM(label),'shape=polygon, sides=3]'
    CASE(status_deleted_result)
      WRITE(file_ptr,*) TRIM(label),'style=filled, shape=polygon, sides=5]'
    CASE(status_calc_approx_ener)
      WRITE(file_ptr,*) TRIM(label),'style=filled, color=brown]'
    CASE(status_calculate_energy)
      WRITE(file_ptr,*) TRIM(label),'style=filled, color=wheat]'
    CASE(status_calculate_MD)
      WRITE(file_ptr,*) TRIM(label),'style=filled, color=goldenrod]'
    CASE(status_calculate_NMC_steps)
      WRITE(file_ptr,*) TRIM(label),'style=filled, color=peru]'
    CASE DEFAULT
      CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
              routineP,"element status"//cp_to_string(stat),error,failure)
    END SELECT
    CALL close_file(unit_number=file_ptr, keep_preconnection=.TRUE.)
  END SUBROUTINE write_color

! *****************************************************************************
!> \brief creates an new branch (hence a new element is created)
!> \param parent_nr tree element number of element one level up
!> \param child_nr tree element number of actual element
!> \param acc flag for accepted or not accepted branch (left,right)
!> \param tmc_params param environment for creating the file name
!> \param tree index of the tree (0=global tree)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE create_dot_branch(parent_nr,child_nr,acc,tmc_params,tree,error)
    INTEGER                                  :: parent_nr, child_nr
    LOGICAL                                  :: acc
    TYPE(tmc_param_type), POINTER            :: tmc_params
    INTEGER                                  :: tree
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_dot_branch', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: file_ptr
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

    CALL open_file(file_name=get_dot_file_name(tmc_params, tree, error), &
                   file_status="OLD", file_action="WRITE", &
                   file_position="APPEND", unit_number=file_ptr)
    IF(acc) THEN
      WRITE(file_ptr,*)parent_nr," -> ", child_nr, ":nw [color=darkolivegreen1]"
    ELSE
      WRITE(file_ptr,*)parent_nr," -> ", child_nr, ":ne [color=coral]"
    END IF
    CALL close_file(unit_number=file_ptr, keep_preconnection=.TRUE.)
  END SUBROUTINE create_dot_branch

! *****************************************************************************
!> \brief interfaces the creating of a branch for subtree elements
!> \param new_element the actual subtree element
!> \param conf the subtree index and hence the index for filename
!> \param tmc_params ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE create_dot(new_element, conf, tmc_params, error)
    TYPE(tree_type), POINTER                 :: new_element
    INTEGER                                  :: conf
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_dot', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(new_element),cp_failure_level,routineP,error,failure)
    CPPrecondition(conf.GT.0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

    CALL create_dot_branch(parent_nr=new_element%parent%nr,&
                           child_nr=new_element%nr,&
                           acc=ASSOCIATED(new_element%parent%acc,new_element), &
                           tmc_params=tmc_params, tree=conf, error=error)
  END SUBROUTINE create_dot

! *****************************************************************************
!> \brief creates new dot and arrow from element one level up (for subtree)
!>        additional handling of nodes with swaped elements
!> \param new_element the actual global element
!> \param tmc_params ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE create_global_tree_dot(new_element, tmc_params, error)
    TYPE(global_tree_type), POINTER          :: new_element
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_global_tree_dot', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=1000)                      :: list_of_nr
    INTEGER                                  :: file_ptr, i, ref_count
    LOGICAL                                  :: failure
    TYPE(gt_elem_list_type), POINTER         :: tmp_pt_list_elem

    failure = .FALSE.
    NULLIFY(tmp_pt_list_elem)

    CPPrecondition(ASSOCIATED(new_element),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

    ! creating list with configuration numbers (of subtrees)
    list_of_nr = ""
    ! the order of subtrees
    DO i=1,SIZE(new_element%conf(:))
      WRITE(list_of_nr, *)TRIM(ADJUSTL(list_of_nr)),new_element%conf(i)%elem%sub_tree_nr
    END DO
    ! the used subtree elements
    WRITE(list_of_nr, *)TRIM(ADJUSTL(list_of_nr)),'\n '
    DO i=1,SIZE(new_element%conf(:))
      WRITE(list_of_nr, *)TRIM(ADJUSTL(list_of_nr))," ", new_element%conf(i)%elem%nr
    END DO
    ! print out the references of each subtree element
    IF(DEBUG.GT.8) THEN
      WRITE(list_of_nr, *)TRIM(ADJUSTL(list_of_nr)),'\n ref'
      DO i=1,SIZE(new_element%conf(:))
        ref_count = 0
        tmp_pt_list_elem => new_element%conf(i)%elem%gt_nodes_references
        DO WHILE(ASSOCIATED(tmp_pt_list_elem))
           ref_count = ref_count +1
           ! create a list with all references
           IF(.FALSE.) WRITE(list_of_nr, *)TRIM(ADJUSTL(list_of_nr))," ", tmp_pt_list_elem%gt_elem%nr
           tmp_pt_list_elem => tmp_pt_list_elem%next
        END DO
        ! print a list with all references
        IF(.FALSE.) WRITE(list_of_nr, *)TRIM(ADJUSTL(list_of_nr)),' | '
        ! print only the amount of references
        IF(.TRUE.)  WRITE(list_of_nr, *) TRIM(ADJUSTL(list_of_nr)), ref_count,' | '
      END DO
    END IF

    IF(.NOT.ASSOCIATED(new_element%parent)) THEN
      CALL cp_assert(new_element%nr.LE.1,&
                    cp_warning_level,cp_assertion_failed,routineP,&
                    "try to create dot, but no parent on node "//&
                    cp_to_string(new_element%nr)//"exists",&
                    error)
    ELSE
      CALL create_dot_branch(parent_nr=new_element%parent%nr,&
                             child_nr=new_element%nr,&
                             acc=ASSOCIATED(new_element%parent%acc,new_element),&
                             tmc_params=tmc_params, tree=0, error=error)
    END IF
    ! write in dot file
    CALL open_file(file_name=get_dot_file_name(tmc_params, 0, error), &
                   file_status="OLD", file_action="WRITE", &
                   file_position="APPEND", unit_number=file_ptr)
    IF(new_element%swaped)THEN
      WRITE(file_ptr,*) new_element%nr,'[label="',new_element%nr, ' |', new_element%mv_conf,' |', &
                        mv_type_swap_conf,'\n ',&
                        TRIM(ADJUSTL(list_of_nr)), '", shape=polygon, peripheries=3, sides=5]'
    ELSE
      WRITE(file_ptr,*) new_element%nr,'[label="',new_element%nr, ' |', new_element%mv_conf,' |', &
                        new_element%conf(new_element%mv_conf)%elem%move_type,'\n ',&
                        TRIM(ADJUSTL(list_of_nr)), '"]'
    END IF
    CALL close_file(file_ptr, keep_preconnection=.TRUE.)
  END SUBROUTINE create_global_tree_dot

! *****************************************************************************
!> \brief interfaces the change of color for subtree elements
!>        on the basis of the element status
!> \param tree_element the actual global element
!> \param tmc_params ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE create_dot_color(tree_element, tmc_params, error)
    TYPE(tree_type), POINTER                 :: tree_element
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_dot_color', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=1000)                      :: list_of_nr
    INTEGER                                  :: ref_count
    LOGICAL                                  :: failure
    TYPE(gt_elem_list_type), POINTER         :: tmp_pt_list_elem

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(tree_element),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

    IF(DEBUG.GT.8)THEN
      list_of_nr = ""
      tmp_pt_list_elem => tree_element%gt_nodes_references
      ref_count = 0
      DO WHILE(ASSOCIATED(tmp_pt_list_elem))
         ref_count = ref_count +1
         ! print a list with all references
         IF(.FALSE.) THEN
            WRITE(list_of_nr, *)TRIM(ADJUSTL(list_of_nr))," ", tmp_pt_list_elem%gt_elem%nr
            WRITE(list_of_nr, *)TRIM(ADJUSTL(list_of_nr)),' | '
         END IF
         ! print only the amount of references
         IF(.TRUE.)  WRITE(list_of_nr, *) ref_count,' | '
         tmp_pt_list_elem => tmp_pt_list_elem%next
      END DO
      WRITE(*,*)"mark subtree",tree_element%sub_tree_nr," node", tree_element%nr, " with status ", &
            tree_element%stat, "ref ", TRIM(ADJUSTL(list_of_nr))
    END IF

    CALL write_color(node_nr=tree_element%nr, stat=tree_element%stat, &
             filename=get_dot_file_name(tmc_params,tree_element%sub_tree_nr,error), error=error)
  END SUBROUTINE create_dot_color

! *****************************************************************************
!> \brief interfaces the change of color for global tree  node
!>        on the basis of the element status
!> \param gt_tree_element the actual global element
!> \param tmc_params ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE create_global_tree_dot_color(gt_tree_element, tmc_params, error)
    TYPE(global_tree_type), POINTER          :: gt_tree_element
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_global_tree_dot_color', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(gt_tree_element),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

     IF(DEBUG.GT.8) WRITE(*,*)"mark global tree node color", gt_tree_element%nr, gt_tree_element%stat
     CALL write_color(node_nr=gt_tree_element%nr, stat=gt_tree_element%stat, &
                      filename=get_dot_file_name(tmc_params,0,error), error=error)
   END SUBROUTINE create_global_tree_dot_color

!! *****************************************************************************
!!> \brief prints out dot file for a whole subtree below the entered element
!!> \param current the actual subtree element
!!> \param conf index of the subtree
!!> \param error variable to control error logging, stopping,...
!!>        see module cp_error_handling
!!> \author Mandes 12.2012
!! *****************************************************************************
!  RECURSIVE SUBROUTINE create_tree(current, conf, filename, error)
!    TYPE (tree_type), POINTER                :: current
!    INTEGER                                  :: conf
!    CHARACTER(LEN=*)                         :: filename
!    TYPE(cp_error_type), INTENT(inout)       :: error
!
!    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_tree', &
!      routineP = moduleN//':'//routineN
!
!    CALL create_dot_color(current, tmc_params, error)
!    IF(ASSOCIATED(current%acc))THEN
!       CALL create_dot_branch(parent_nr=current%nr, child_nr=current%acc%nr, &
!                              acc=.TRUE.,tmc_params=tmc_params, file_single_tree_ptr)
!       WRITE(file_single_tree_ptr,*)current%nr,'[label="', current%nr,"\n ",&
!                                    current%pos(1),"\n ", current%potential,'"]'
!       CALL create_tree(current%acc, conf)
!    END IF
!    IF(ASSOCIATED(current%nacc))THEN
!       CALL create_dot_branch(current%nr,current%acc%nr,.FALSE.,file_single_tree_ptr)
!       WRITE(file_single_tree_ptr,*)current%nr,'[label="', current%nr,"\n ",&
!                                    current%pos(1),"\n ", current%potential,'"]'
!       CALL create_tree(current%nacc, conf)
!    END IF
!  END SUBROUTINE create_tree
END MODULE tmc_dot_tree
