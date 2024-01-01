!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
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
      INTEGER                                          :: ref_count = -1
      TYPE(routine_map_type)                           :: routine_names = routine_map_type()
      TYPE(list_routinestat_type)                      :: routine_stats = list_routinestat_type()
      TYPE(list_callstackentry_type)                   :: callstack = list_callstackentry_type()
      TYPE(callgraph_type)                             :: callgraph = callgraph_type()
      INTEGER                                          :: trace_max = -1
      INTEGER                                          :: trace_unit = -1
      CHARACTER(len=13)                                :: trace_str = ""
      LOGICAL                                          :: trace_all = .FALSE.
   END TYPE timer_env_type

   PUBLIC :: timer_env_type

END MODULE timings_types

! **************************************************************************************************

