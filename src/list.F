!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!


! *****************************************************************************
!> \brief A common interface for all instances of list__valuetype_.template .
!> \par History
!>      12.2012 first version [ole]
!> \author Ole Schuett
! *****************************************************************************
	
MODULE list
  USE list_callstackentry
  USE list_routinestat
  USE list_timerenv

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'list'
  
  PUBLIC :: list_init, list_isready, list_push, list_pop, list_peek,&
   list_insert, list_set, list_get, list_del, list_clear, list_size, list_destroy



  INTERFACE list_init
     MODULE PROCEDURE list_timerenv_init, list_routinestat_init, &
       list_callstackentry_init
  END INTERFACE
  
  INTERFACE list_isready
     MODULE PROCEDURE list_timerenv_isready, &
        list_routinestat_isready, list_callstackentry_isready
  END INTERFACE
  
  INTERFACE list_push
     MODULE PROCEDURE list_timerenv_push, list_routinestat_push, &
       list_callstackentry_push
  END INTERFACE
  
  INTERFACE list_pop
     MODULE PROCEDURE list_timerenv_pop, list_routinestat_pop, &
       list_callstackentry_pop
  END INTERFACE
  
  INTERFACE list_peek
     MODULE PROCEDURE list_timerenv_peek, list_routinestat_peek, &
       list_callstackentry_peek
  END INTERFACE
    
  INTERFACE list_insert
     MODULE PROCEDURE list_timerenv_insert, & 
        list_routinestat_insert, list_callstackentry_insert
  END INTERFACE

  INTERFACE list_set
     MODULE PROCEDURE list_timerenv_set, list_routinestat_set, &
       list_callstackentry_set
  END INTERFACE
  
  INTERFACE list_get
     MODULE PROCEDURE list_timerenv_get, list_routinestat_get, &
       list_callstackentry_get
  END INTERFACE
  
  INTERFACE list_del
     MODULE PROCEDURE list_timerenv_del, list_routinestat_del, &
       list_callstackentry_del
  END INTERFACE
  
  INTERFACE list_clear
     MODULE PROCEDURE list_timerenv_clear, list_routinestat_clear,&
       list_callstackentry_clear
  END INTERFACE
  
  INTERFACE list_size
     MODULE PROCEDURE list_timerenv_size, list_routinestat_size, &
       list_callstackentry_size
  END INTERFACE
  
  INTERFACE list_destroy
     MODULE PROCEDURE list_timerenv_destroy, &
        list_routinestat_destroy, list_callstackentry_destroy
  END INTERFACE
    
END MODULE list
