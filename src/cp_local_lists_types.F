!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_local_lists_types [1.0] *
!!
!!   NAME
!!     cp_local_lists_types
!!
!!   FUNCTION
!!     stores a lists of integer that are local to a processor.
!!     The idea is that these integers represent objects that are distributed
!!     between the different processors.
!!     The ordering is just to make some operation more efficent, logically
!!     these lists are like sets.
!!     Some operations assume that the integers represent a range of values
!!     from 1 to a (not too big) maxval, and that an element is present just
!!     once, and only on a processor (these assumption are marked in the 
!!     documentation of such operations).
!!     The concrete task for which this structure was developed was 
!!     distributing atoms between the processors.
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE cp_local_lists_types
  USE cp_array_i_utils,                ONLY: cp_1d_i_bsearch,&
                                             cp_1d_i_p_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_check,&
                                             cp_error_dealloc_ref,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: wp=>dp
  USE message_passing,                 ONLY: mp_max,&
                                             mp_sum
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_local_lists_types'
  INTEGER, SAVE, PRIVATE :: cp_ll_last_id_nr=0

  PUBLIC :: cp_local_lists_type, cp_local_lists_p_type
  PUBLIC :: cp_ll_create, cp_ll_retain, cp_ll_release, &
       cp_ll_build_global_dist, cp_ll_destroy_global_dist,&
       cp_ll_search_el, cp_ll_sort_list, cp_ll_write
!***
!****************************************************************************

!!****s* cp_local_lists_types/cp_local_lists_type [1.0] *
!!
!!   NAME
!!     cp_local_lists_type
!!
!!   FUNCTION
!!     structure to store local (to a processor) ordered lists of integers.
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - ref_count: reference count (see doc/ReferenceCounting.html)
!!     - n_el: n_el(i) is number of elements of list(i)
!!     - list: list(i) contains an ordered list of integer (the array
!!       might be bigger than n_el(i), but the extra elements should be
!!       ignored)
!!     - para_env: the parallel environment underlying the distribution
!!     - global_distribution: global_distribution(i) is the id of the
!!       processor that has "i". This is allocated only with
!!       cp_ll_build_global_distribution, try to avoid it (uses more
!!       memory)
!!     - listbased_distribution: true if each list has its own 
!!       distribution
!!       
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     06.2002 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE cp_local_lists_type
     INTEGER :: ref_count, id_nr
     LOGICAL :: listbased_distribution
     INTEGER, DIMENSION(:), POINTER :: n_el
     TYPE(cp_1d_i_p_type), DIMENSION(:), POINTER :: list, global_distribution
     TYPE(cp_para_env_type), POINTER :: para_env
  END TYPE cp_local_lists_type
!!***
  !****************************************************************************

!!****s* cp_local_lists_types/cp_local_lists_p_type [1.0] *
!!
!!   NAME
!!     cp_local_lists_p_type
!!
!!   FUNCTION
!!     pointer to a cp_local_lists_type
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - local_lists: the pointer to the local_lists
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     06.2003 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE cp_local_lists_p_type
     TYPE(cp_local_lists_type), POINTER :: local_lists
  END TYPE cp_local_lists_p_type
!!***
  !****************************************************************************

CONTAINS

!!****f* cp_local_lists_types/cp_ll_create [1.0] *
!!
!!   NAME
!!     cp_ll_create
!!
!!   SYNOPSIS
!!     Subroutine cp_ll_create(local_lists, para_env,&
!!         listbased_distribution, n_el, n_lists, error)
!!       Type(cp_local_lists_type), Pointer:: local_lists
!!       Integer, Intent (IN), Optional:: n_lists
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Logical, Intent (IN), Optional:: listbased_distribution
!!       Integer, Dimension(:), Intent (IN), Optional:: n_el
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_ll_create
!!
!!   FUNCTION
!!     creates a local list
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - local_lists: the lists to create
!!     - n_lists: number of lists to create (defaults to 1, or size(n_el))
!!     - para_env: parallel environement to be used
!!     - listbased_distribution: if each list has its own distribution
!!       (defaults to false)
!!     - n_el: number of elements in each list (defaults to 0)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_ll_create(local_lists,para_env,listbased_distribution,&
     n_el,n_lists,error)
  TYPE(cp_local_lists_type), POINTER :: local_lists
  INTEGER, INTENT(in), OPTIONAL :: n_lists
  TYPE(cp_para_env_type), POINTER :: para_env
  LOGICAL, INTENT(in), OPTIONAL :: listbased_distribution
  INTEGER, DIMENSION(:), INTENT(in), OPTIONAL :: n_el
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_ll_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat, ilist, my_n_lists

  failure=.FALSE.
  my_n_lists=1
  IF(PRESENT(n_el)) my_n_lists=SIZE(n_el)
  IF(PRESENT(n_lists)) my_n_lists=n_lists
  
  ALLOCATE(local_lists,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     cp_ll_last_id_nr=cp_ll_last_id_nr+1
     local_lists%id_nr=cp_ll_last_id_nr

     local_lists%ref_count=1
     
     local_lists%para_env => para_env
     CALL cp_para_env_retain(para_env,error=error)
     
     NULLIFY(local_lists%global_distribution)
     
     local_lists%listbased_distribution=.FALSE.
     IF (PRESENT(listbased_distribution)) &
          local_lists%listbased_distribution=listbased_distribution
     
     ALLOCATE(local_lists%n_el(my_n_lists), local_lists%list(my_n_lists),&
          stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

     IF (PRESENT(n_el)) THEN
        local_lists%n_el=n_el(1:my_n_lists)
     ELSE
        local_lists%n_el=0
     END IF

     DO ilist=1,my_n_lists
        ALLOCATE(local_lists%list(ilist)%array(n_el(ilist)),stat=stat)
        local_lists%list(ilist)%array=-1
     END DO

  END IF
END SUBROUTINE cp_ll_create
!***************************************************************************

!!****f* cp_local_lists_types/cp_ll_retain [1.0] *
!!
!!   NAME
!!     cp_ll_retain
!!
!!   SYNOPSIS
!!     Subroutine cp_ll_retain(local_lists, error)
!!       Type(cp_local_lists_type), Pointer:: local_lists
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_ll_retain
!!
!!   FUNCTION
!!     retains a local_lists
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - local_lists: the object to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_ll_retain(local_lists,error)
  TYPE(cp_local_lists_type), POINTER :: local_lists
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_ll_retain',&
        routineP=moduleN//':'//routineN

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(local_lists),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(local_lists%ref_count>0,cp_failure_level,routineP,error)
     local_lists%ref_count=local_lists%ref_count+1
  END IF
END SUBROUTINE cp_ll_retain
!***************************************************************************

!!****f* cp_local_lists_types/cp_ll_release [1.0] *
!!
!!   NAME
!!     cp_ll_release
!!
!!   SYNOPSIS
!!     Subroutine cp_ll_release(local_lists, error)
!!       Type(cp_local_lists_type), Pointer:: local_lists
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_ll_release
!!
!!   FUNCTION
!!     releases the given local_lists
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - local_lists: the object to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_ll_release(local_lists,error)
  TYPE(cp_local_lists_type), POINTER :: local_lists
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_ll_release',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat, ilist

  failure=.FALSE.
  
  IF (ASSOCIATED(local_lists)) THEN
     CPPreconditionNoFail(local_lists%ref_count>0,cp_failure_level,routineP,error)
     local_lists%ref_count=local_lists%ref_count-1
     IF (local_lists%ref_count==0) THEN
        local_lists%ref_count=1
        CALL cp_ll_destroy_global_dist(local_lists,error=error)
        local_lists%ref_count=0

        DEALLOCATE(local_lists%n_el,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        
        DO ilist=1,SIZE(local_lists%list)
           DEALLOCATE(local_lists%list(ilist)%array,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END DO
        DEALLOCATE(local_lists%list,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        
        CALL cp_para_env_release(local_lists%para_env,error=error)

        DEALLOCATE(local_lists,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
END SUBROUTINE cp_ll_release
!***************************************************************************

!!****f* cp_local_lists_types/cp_ll_sort_list [1.0] *
!!
!!   NAME
!!     cp_ll_sort_list
!!
!!   SYNOPSIS
!!     Subroutine cp_ll_sort_list(local_lists, i_list, error)
!!       Type(cp_local_lists_type), Pointer:: local_lists
!!       Integer, Intent (IN), Optional:: i_list
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_ll_sort_list
!!
!!   FUNCTION
!!     sorts the given list (if not specified all)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - local_lists: the lists to sort
!!     - i_list: the index of the list to sort (if not specified sorts
!!       all the lists)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_ll_sort_list(local_lists,i_list,error)
  TYPE(cp_local_lists_type), POINTER :: local_lists
  INTEGER, INTENT(in), OPTIONAL :: i_list
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_ll_sort_list',&
        routineP=moduleN//':'//routineN
  INTEGER, DIMENSION(:), ALLOCATABLE :: indexes
  INTEGER :: stat,ilist

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(local_lists),cp_failure_level,routineP,error,failure)
  CPPrecondition(local_lists%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (PRESENT(i_list)) THEN
        ALLOCATE(indexes(local_lists%n_el(i_list)),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

        CALL sort(local_lists%list(i_list)%array,&
             local_lists%n_el(i_list),indexes)

        DEALLOCATE(indexes,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     ELSE
        DO ilist=1,SIZE(local_lists%n_el)
           ALLOCATE(indexes(local_lists%n_el(ilist)),stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           
           CALL sort(local_lists%list(ilist)%array,&
                local_lists%n_el(ilist),indexes)
           
           DEALLOCATE(indexes,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END DO
     END IF
  END IF
END SUBROUTINE cp_ll_sort_list
!***************************************************************************

!!****f* cp_local_lists_types/cp_ll_search_el [1.0] *
!!
!!   NAME
!!     cp_ll_search_el
!!
!!   SYNOPSIS
!!     Subroutine cp_ll_search_el(local_lists, element, found,&
!!         i_list_to_search, i_pos, i_list, error)
!!       Type(cp_local_lists_type), Pointer:: local_lists
!!       Integer, Intent (IN):: element
!!       Logical, Intent (OUT):: found
!!       Integer, Intent (IN), Optional:: i_list_to_search
!!       Integer, Optional, Intent (OUT):: i_pos, i_list
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_ll_search_el
!!
!!   FUNCTION
!!     returns true if the given element is contained in the requested list
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - local_lists: the lists to be looked in
!!     - element: the element to look for
!!     - found: returns true if the element was found
!!     - i_list_to_search: index of the list to be searched (if not given
!!       all the lists are searched)
!!     - i_pos: optional, returns the position of the element in the list
!!       (makes sense only if an element was found, or if you searched only
!!       a single list, in that case it is the place where you should insert
!!       the element)
!!     - i_list: optional, returns the index of the list where the element
!!       was found (makes sense only if an element was found)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     06.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_ll_search_el(local_lists, element, found, i_list_to_search, &
     i_pos,i_list, error)
  TYPE(cp_local_lists_type), POINTER :: local_lists
  INTEGER, INTENT(in) :: element
  LOGICAL, INTENT(out) :: found
  INTEGER, INTENT(in), OPTIONAL :: i_list_to_search
  INTEGER, OPTIONAL, INTENT(out) :: i_pos, i_list
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_ll_search_el',&
       routineP=moduleN//':'//routineN
  INTEGER :: index, ilist
  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(local_lists),cp_failure_level,routineP,error,failure)
  CPPrecondition(local_lists%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     found=.FALSE.
     IF (PRESENT(i_list_to_search)) THEN
        ilist=i_list_to_search
        index=cp_1d_i_bsearch(array=local_lists%list(ilist)%array,&
             el=element, l_index=1, u_index=local_lists%n_el(ilist),&
             error=error)
        IF (index <= local_lists%n_el(ilist)) THEN
           found=(local_lists%list(ilist)%array(index)==element)
        END IF
     ELSE
        DO ilist=1,SIZE(local_lists%list)
           index=cp_1d_i_bsearch(array=local_lists%list(ilist)%array,&
                el=element, l_index=1, u_index=local_lists%n_el(ilist),&
                error=error)
           IF (index <= local_lists%n_el(ilist)) THEN
              found=(local_lists%list(ilist)%array(index)==element)
              IF (found) EXIT
           END IF
        END DO
     END IF
     IF (PRESENT(i_list)) i_list=ilist
     IF (PRESENT(i_pos)) i_pos=index
  END IF
END SUBROUTINE cp_ll_search_el
!***************************************************************************

!!****f* cp_local_lists_types/cp_ll_build_global_dist [1.0] *
!!
!!   NAME
!!     cp_ll_build_global_dist
!!
!!   SYNOPSIS
!!     Subroutine cp_ll_build_global_dist(local_lists, error)
!!       Type(cp_local_lists_type), Pointer:: local_lists
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_ll_build_global_dist
!!
!!   FUNCTION
!!     build an array with the id of the processor that owns that index
!!     (-1 if no processor owns it)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - local_lists: the list of which you want to build the global
!!       distribution
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     06.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_ll_build_global_dist(local_lists,error)
  TYPE(cp_local_lists_type), POINTER :: local_lists
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_ll_build_global_dist',&
       routineP=moduleN//':'//routineN
  INTEGER :: ilist, stat, max_el, max_loc, max_el1, i, my_id
  INTEGER, DIMENSION(:), POINTER :: g_dist

  failure=.FALSE.
  NULLIFY(g_dist)

  CPPrecondition(ASSOCIATED(local_lists),cp_failure_level,routineP,error,failure)
  CPPrecondition(local_lists%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (local_lists%listbased_distribution) THEN

        IF (.NOT.ASSOCIATED(local_lists%global_distribution)) THEN
           ALLOCATE(local_lists%global_distribution(SIZE(local_lists%list)),&
                stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           DO ilist=1,SIZE(local_lists%list)
              NULLIFY(local_lists%global_distribution(ilist)%array)
           END DO
        END IF

        DO ilist=1,SIZE(local_lists%list)
           max_el=0
           !$omp parallel do private(i,max_loc) reduction(max:max_el)
           DO i=1,local_lists%n_el(ilist)
              max_loc=local_lists%list(ilist)%array(i)
              IF (max_el < max_loc) max_el=max_loc
           END DO
           CALL mp_max(max_el, local_lists%para_env%group)

           IF (ASSOCIATED(local_lists%global_distribution(ilist)%array)) THEN
              IF (SIZE(local_lists%global_distribution(ilist)%array)<&
                   max_el) THEN
                 DEALLOCATE(local_lists%global_distribution(ilist)%array,&
                      stat=stat)
                 CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
              END IF
           END IF
           IF (.NOT.ASSOCIATED(local_lists%global_distribution &
                (ilist)%array)) THEN
              ALLOCATE(local_lists%global_distribution(ilist)%array(max_el),&
                   stat=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           END IF

           g_dist => local_lists%global_distribution(ilist)%array
           g_dist=0
           my_id=local_lists%para_env%mepos+1

           !$omp parallel do
           DO i=1,local_lists%n_el(ilist)
              g_dist(local_lists%list(ilist)%array(i))=my_id
           END DO
           CALL mp_sum(g_dist,local_lists%para_env%group)
           !$omp parallel do
           DO i=1,SIZE(g_dist)
              g_dist(i)=g_dist(i)-1
           END DO

           my_id=my_id-1
           !$omp parallel do
           DO i=1,local_lists%n_el(ilist)
              IF (g_dist(local_lists%list(ilist)%array(i))/=my_id) THEN
                 CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                      routineP,"invalid distribution, two processors have "//&
                      "the same number"//&
                      CPSourceFileRef,&
                      error,failure)
              END IF
           END DO
        END DO
     ELSE
        IF (.NOT.ASSOCIATED(local_lists%global_distribution)) THEN
           ALLOCATE(local_lists%global_distribution(1),&
                stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           NULLIFY(local_lists%global_distribution(1)%array)
        END IF
        max_el=0
        DO ilist=1,SIZE(local_lists%list)
           max_el1=0
           !$omp parallel do private(i,max_loc) reduction(max:max_el1)
           DO i=1,local_lists%n_el(ilist)
              max_loc=local_lists%list(ilist)%array(i)
              IF (max_el1 < max_loc) max_el1=max_loc
           END DO
           IF (max_el < max_el1) max_el=max_el1
        END DO
        CALL mp_max(max_el, local_lists%para_env%group)
        IF (ASSOCIATED(local_lists%global_distribution(1)%array)) THEN
           IF (SIZE(local_lists%global_distribution(1)%array)<&
                max_el) THEN
              DEALLOCATE(local_lists%global_distribution(1)%array,&
                   stat=stat)
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           END IF
        END IF
        IF (.NOT.ASSOCIATED(local_lists%global_distribution &
             (1)%array)) THEN
           ALLOCATE(local_lists%global_distribution(1)%array(max_el),&
                stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END IF

        g_dist => local_lists%global_distribution(1)%array
        g_dist=0
        my_id=local_lists%para_env%mepos+1

        DO ilist=1,SIZE(local_lists%list)
           !$omp parallel do
           DO i=1,local_lists%n_el(ilist)
              ! use += to catch multiple int on same proc?
              g_dist(local_lists%list(ilist)%array(i))=my_id
           END DO
        END DO
        CALL mp_sum(g_dist,local_lists%para_env%group)
        !$omp parallel do
        DO i=1,SIZE(g_dist)
           g_dist(i)=g_dist(i)-1
        END DO

        my_id=my_id-1
        DO ilist=1,SIZE(local_lists%list)
           !$omp parallel do
           DO i=1,local_lists%n_el(ilist)
              IF (g_dist(local_lists%list(ilist)%array(i))/=my_id) THEN
                 CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                      routineP,"invalid distribution, two processors have "//&
                      "the same number (list="//cp_to_string(ilist)//", i="//&
                      cp_to_string(i)//",g_index="//&
                      cp_to_string(local_lists%list(ilist)%array(i))//")"//&
                      CPSourceFileRef,&
                      error,failure)
              END IF
           END DO
        END DO
     END IF
  END IF
END SUBROUTINE cp_ll_build_global_dist
!***************************************************************************

!!****f* cp_local_lists_types/cp_ll_destroy_global_dist [1.0] *
!!
!!   NAME
!!     cp_ll_destroy_global_dist
!!
!!   SYNOPSIS
!!     Subroutine cp_ll_destroy_global_dist(local_lists, error)
!!       Type(cp_local_lists_type), Pointer:: local_lists
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_ll_destroy_global_dist
!!
!!   FUNCTION
!!     deallocates the global dist array
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - local_lists: the list of which you want to destroy the global
!!       distribution
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     06.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_ll_destroy_global_dist(local_lists,error)
  TYPE(cp_local_lists_type), POINTER :: local_lists
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_ll_destroy_global_dist',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat, ilist

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(local_lists),cp_failure_level,routineP,error,failure)
  CPPrecondition(local_lists%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (ASSOCIATED(local_lists%global_distribution)) THEN
        DO ilist=1,SIZE(local_lists%global_distribution)
           IF (ASSOCIATED(local_lists%global_distribution(ilist)%array)) THEN
              DEALLOCATE(local_lists%global_distribution(ilist)%array,&
                   stat=stat)
           END IF
        END DO

        DEALLOCATE(local_lists%global_distribution, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
END SUBROUTINE cp_ll_destroy_global_dist
!***************************************************************************

!!****f* cp_local_lists_type/cp_ll_write [1.0] *
!!
!!   NAME
!!     cp_ll_write
!!
!!   FUNCTION
!!     writes a description of the local distribution to the given unit
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - local_lists: the list to describe
!!     - unit_nr: the unit to write to
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     06.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_ll_write(local_lists, unit_nr, error)
  TYPE(cp_local_lists_type), pointer :: local_lists
  INTEGER, INTENT(in) :: unit_nr
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_ll_write',&
        routineP=moduleN//':'//routineN
  INTEGER :: ilist, i

  failure=.false.

  IF (ASSOCIATED(local_lists)) THEN
     WRITE(unit=unit_nr, fmt="(' <local_lists>{ id_nr=',i8,',ref_count=',i8,',')")&
          local_lists%id_nr, local_lists%ref_count
     IF (ASSOCIATED(local_lists%n_el)) THEN
        WRITE(unit=unit_nr, fmt="(' n_el=(')", advance="no")
        DO i=1,SIZE(local_lists%n_el)
           WRITE(unit=unit_nr, fmt="(i6,',')", advance="no") local_lists%n_el(i)
        END DO
        WRITE(unit=unit_nr, fmt="('),')")
     ELSE
        WRITE(unit=unit_nr, fmt="(' n_el=*null*,')")
     END IF
     IF (ASSOCIATED(local_lists%list)) THEN
        WRITE(unit=unit_nr, fmt="(' list=(')")
        DO ilist=1,size(local_lists%list)
           IF (ASSOCIATED(local_lists%list(ilist)%array)) THEN
              WRITE(unit=unit_nr, fmt="('  (')", advance="no")
              DO i=1,SIZE(local_lists%list(ilist)%array)
                 WRITE(unit=unit_nr, fmt="(i6,',')", advance="no")&
                      local_lists%list(ilist)%array(i)
              END DO
              WRITE(unit=unit_nr, fmt="('),')")
           ELSE
              WRITE(unit=unit_nr, fmt="('  *null*,')")
           END IF
        END DO
        WRITE(unit=unit_nr, fmt="(' )')")
     ELSE
        WRITE(unit=unit_nr, fmt="(' list=*null*,')")
     END IF
     WRITE(unit=unit_nr, fmt="(' }')")
  ELSE
     WRITE(unit=unit_nr, fmt="(' <local_lists *null*>,')")
  END IF
END SUBROUTINE cp_ll_write
!***************************************************************************

END MODULE cp_local_lists_types
