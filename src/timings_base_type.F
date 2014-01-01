!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief Types used by timings.F and timings_report.F
!>        The types in this module are used within dict or list, which are 
!>	  in turn used in timer_env_type from timings_type.F
!>        Due to the fortran restriction on cicular module-dependencies these
!>        "inner-types" had to go into a separate module. 
!> \par History
!>      12.2012 Created [ole]	
!> \author Ole Schuett
! *****************************************************************************
MODULE timings_base_type       
 
  USE kinds,                           ONLY: default_string_length,&
                                             dp

  IMPLICIT NONE
  PRIVATE
  
  TYPE routine_stat_type
     INTEGER       :: routine_id
     CHARACTER(len=default_string_length) :: routineN
     REAL(kind=dp) :: excl_walltime_accu
     REAL(kind=dp) :: incl_walltime_accu
     REAL(kind=dp) :: excl_energy_accu
     REAL(kind=dp) :: incl_energy_accu
     INTEGER       :: active_calls
     INTEGER       :: total_calls
     INTEGER       :: stackdepth_accu
     LOGICAL       :: trace
  END TYPE routine_stat_type
  
  TYPE call_stat_type
     INTEGER       :: total_calls
     REAL(kind=dp) :: incl_walltime_accu
     REAL(kind=dp) :: incl_energy_accu
  END TYPE call_stat_type

  
  TYPE callstack_entry_type
     INTEGER       :: routine_id
     REAL(kind=dp) :: walltime_start
     REAL(kind=dp) :: energy_start
  END TYPE callstack_entry_type
  
  PUBLIC :: routine_stat_type, call_stat_type, callstack_entry_type

END MODULE timings_base_type   

! *****************************************************************************

