!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types used by timings.F and timings_report.F
!>        Due to the fortran restriction on cicular module-dependencies 
!>        the types, which are used through dict or list had to go 
!>        into the separate module timings_base_type.F 
!> \par History
!>      12.2012 Created [ole]	
!> \author Ole Schuett
! *****************************************************************************
MODULE timings_types      
  USE dict_i4tuple_callstat,           ONLY: dict_i4tuple_callstat_type
  USE dict_str_i4,                     ONLY: dict_str_i4_type
  USE list_callstackentry,             ONLY: list_callstackentry_type
  USE list_routinestat,                ONLY: list_routinestat_type

  IMPLICIT NONE
  PRIVATE
  
  TYPE timer_env_type
     INTEGER                                          :: ref_count
     TYPE(dict_str_i4_type)                           :: routine_names
     TYPE(list_routinestat_type)                      :: routine_stats
     TYPE(list_callstackentry_type)                   :: callstack
     TYPE(dict_i4tuple_callstat_type)                 :: callgraph
     INTEGER                                          :: trace_max
     INTEGER                                          :: trace_unit
     CHARACTER(len=13)                                :: trace_str
     LOGICAL                                          :: trace_all
  END TYPE timer_env_type                                
 
  PUBLIC :: timer_env_type

END MODULE timings_types   

! *****************************************************************************

