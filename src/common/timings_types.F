!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types used by timings.F and timings_report.F
!>        Due to the fortran restriction on cicular module-dependencies
!>        the types, which are used through dict or list had to go
!>        into the separate module timings_base_type.F
!> \par History
!>      12.2012 Created [ole]
!> \author Ole Schuett
! **************************************************************************************************
MODULE timings_types
   USE callgraph,                       ONLY: callgraph_type
   USE list_callstackentry,             ONLY: list_callstackentry_type
   USE list_routinestat,                ONLY: list_routinestat_type
   USE routine_map,                     ONLY: routine_map_type

   IMPLICIT NONE
   PRIVATE

   TYPE timer_env_type
      INTEGER                                          :: ref_count
      TYPE(routine_map_type)                           :: routine_names
      TYPE(list_routinestat_type)                      :: routine_stats
      TYPE(list_callstackentry_type)                   :: callstack
      TYPE(callgraph_type)                             :: callgraph
      INTEGER                                          :: trace_max
      INTEGER                                          :: trace_unit
      CHARACTER(len=13)                                :: trace_str
      LOGICAL                                          :: trace_all
   END TYPE timer_env_type

   PUBLIC :: timer_env_type

END MODULE timings_types

! **************************************************************************************************

