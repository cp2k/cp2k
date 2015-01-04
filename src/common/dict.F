!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!


! *****************************************************************************
!> \brief A common interface for all instances of dict__keytype___valuetype_.template .
!> \par History
!>      12.2012 first version [ole]
!> \author Ole Schuett
! *****************************************************************************

MODULE dict
  USE dict_i4tuple_callstat,           ONLY: &
       dict_i4tuple_callstat_del, dict_i4tuple_callstat_destroy, &
       dict_i4tuple_callstat_get, dict_i4tuple_callstat_haskey, &
       dict_i4tuple_callstat_init, dict_i4tuple_callstat_isready, &
       dict_i4tuple_callstat_items, dict_i4tuple_callstat_set, &
       dict_i4tuple_callstat_size, dict_i4tuple_callstat_update
  USE dict_str_i4,                     ONLY: &
       dict_str_i4_del, dict_str_i4_destroy, dict_str_i4_get, &
       dict_str_i4_haskey, dict_str_i4_init, dict_str_i4_isready, &
       dict_str_i4_items, dict_str_i4_set, dict_str_i4_size, &
       dict_str_i4_update

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dict'
  
  PUBLIC :: dict_init, dict_items, dict_haskey,&
            dict_set, dict_get, dict_size, dict_destroy

  INTERFACE dict_init
     MODULE PROCEDURE dict_str_i4_init, dict_i4tuple_callstat_init
  END INTERFACE
  
  INTERFACE dict_isready
     MODULE PROCEDURE dict_str_i4_isready, dict_i4tuple_callstat_isready
  END INTERFACE
  
  INTERFACE dict_haskey
     MODULE PROCEDURE dict_str_i4_haskey, dict_i4tuple_callstat_haskey
  END INTERFACE
    
  INTERFACE dict_set
     MODULE PROCEDURE dict_str_i4_set, dict_i4tuple_callstat_set
  END INTERFACE
  
  INTERFACE dict_get
     MODULE PROCEDURE dict_str_i4_get, dict_i4tuple_callstat_get 
  END INTERFACE
  
  INTERFACE dict_del
     MODULE PROCEDURE dict_str_i4_del, dict_i4tuple_callstat_del
  END INTERFACE
  
  INTERFACE dict_items
     MODULE PROCEDURE dict_str_i4_items, dict_i4tuple_callstat_items
  END INTERFACE
    
  INTERFACE dict_size
     MODULE PROCEDURE dict_str_i4_size, dict_i4tuple_callstat_size
  END INTERFACE
  
  INTERFACE dict_update
     MODULE PROCEDURE dict_str_i4_update, dict_i4tuple_callstat_update
  END INTERFACE
  
  INTERFACE dict_destroy
     MODULE PROCEDURE dict_str_i4_destroy, dict_i4tuple_callstat_destroy
  END INTERFACE
    
END MODULE dict
