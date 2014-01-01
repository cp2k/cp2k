!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Timing routines for accounting
!> \par History
!>      02.2004 made a stacked version (of stacks...) [Joost VandeVondele]
!>      11.2004 storable timer_envs (for f77 interface) [fawzi]
!>      10.2005 binary search to speed up lookup in timeset [fawzi]
!>      12.2012 Complete rewrite based on dictionaries. [ole]	
!> \author JGH
! *****************************************************************************
MODULE timings_report
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dict,                            ONLY: dict_get,&
                                             dict_haskey,&
                                             dict_items
  USE dict_i4tuple_callstat,           ONLY: dict_i4tuple_callstat_item_type
  USE f77_blas
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_8
  USE list,                            ONLY: list_get,&
                                             list_size
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum
  USE timings,                         ONLY: get_timer_env
  USE timings_base_type,               ONLY: call_stat_type,&
                                             routine_stat_type
  USE timings_types,                   ONLY: timer_env_type
  USE util,                            ONLY: sort

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: accu_type_time=17, accu_type_energy=18

  PUBLIC :: timings_report_print, timings_report_callgraph, accu_type_time, accu_type_energy

  CONTAINS

! *****************************************************************************
!> \brief Print accumulated information on timers
!> \param para_env is needed to collect statistics from other nodes.
!> \par History
!>      none
!> \author JGH
! *****************************************************************************
  SUBROUTINE timings_report_print(iw,r_timings,sort_by_self_time,accu_type,para_env)
    INTEGER, INTENT(IN)                      :: iw
    REAL(KIND=dp), INTENT(IN)                :: r_timings
    LOGICAL, INTENT(IN)                      :: sort_by_self_time
    INTEGER, INTENT(IN)                      :: accu_type
    TYPE(cp_para_env_type), INTENT(IN)       :: para_env

    CHARACTER(LEN=default_string_length)     :: fmt, routineN
    INTEGER                                  :: decimals, i, j, num_routines, &
                                                routine_id
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: indices
    REAL(KIND=dp)                            :: asd, maxtime, mintime
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:), TARGET                   :: max_ewalltime, max_iwalltime, &
                                                sum_ewalltime, sum_iwalltime
    REAL(KIND=dp), DIMENSION(:), POINTER     :: max_walltime
    TYPE(routine_stat_type), POINTER         :: r_stat
    TYPE(timer_env_type), POINTER            :: timer_env

    NULLIFY (max_walltime)

    timer_env => get_timer_env()
    ! find the number of timing events on the parent
    ! note: we assume that timing events on the parent
    !       include all possible timing events.. i.e. there
    !       are no other timing events on the rest of the group.
    !       If this is the case they will be ignored in the
    !       timing report.

    num_routines = list_size(timer_env%routine_stats)
    CALL mp_bcast(num_routines, para_env%source, para_env%group)

    ALLOCATE(max_iwalltime(num_routines))
    ALLOCATE(sum_iwalltime(num_routines))
    ALLOCATE(max_ewalltime(num_routines))
    ALLOCATE(sum_ewalltime(num_routines))

    ! for each clock_name find the needed quantities over all CPUs.
    DO i=1, num_routines
       IF(para_env%mepos == para_env%source) THEN
                  r_stat => list_get(timer_env%routine_stats, i)
                  routineN = r_stat%routineN
       END IF

       CALL mp_bcast(routineN, para_env%source, para_env%group)
       IF(dict_haskey(timer_env%routine_names, routineN)) THEN
           ! assinged routine_id might vary between nodes 
           routine_id = dict_get(timer_env%routine_names, routineN)
           r_stat => list_get(timer_env%routine_stats, routine_id)
           SELECT CASE(accu_type)
           CASE(accu_type_energy)
               max_iwalltime(i) = r_stat%incl_energy_accu
               sum_iwalltime(i) = r_stat%incl_energy_accu
               max_ewalltime(i) = r_stat%excl_energy_accu
               sum_ewalltime(i) = r_stat%excl_energy_accu
           CASE(accu_type_time)
               max_iwalltime(i) = r_stat%incl_walltime_accu
               sum_iwalltime(i) = r_stat%incl_walltime_accu
               max_ewalltime(i) = r_stat%excl_walltime_accu
               sum_ewalltime(i) = r_stat%excl_walltime_accu
           CASE DEFAULT
               STOP "BUG"
           END SELECT
       ELSE
           ! this node did not register the routine 
           max_iwalltime(i) = 0.0_dp
           sum_iwalltime(i) = 0.0_dp
           max_ewalltime(i) = 0.0_dp
           sum_ewalltime(i) = 0.0_dp
       END IF

       CALL mp_max(max_iwalltime(i), para_env%group)
       CALL mp_sum(sum_iwalltime(i), para_env%group)
       CALL mp_max(max_ewalltime(i), para_env%group)
       CALL mp_sum(sum_ewalltime(i), para_env%group)
    ENDDO


    IF(num_routines>0 .AND. iw>0) THEN
       IF (sort_by_self_time) THEN
          max_walltime => max_ewalltime
       ELSE
          max_walltime => max_iwalltime
       END IF

       ALLOCATE(indices(num_routines))
       CALL sort(max_walltime, num_routines, indices)
       maxtime = max_walltime(num_routines)
       mintime = maxtime*r_timings

       WRITE (UNIT=iw,FMT="(/,T2,A)") REPEAT("-",79)
       WRITE (UNIT=iw,FMT="(T2,A,T80,A)") "-","-"
       SELECT CASE(accu_type)
       CASE(accu_type_energy)
           WRITE (UNIT=iw,FMT="(T2,A,T35,A,T80,A)") "-","E N E R G Y","-"
       CASE(accu_type_time)
           WRITE (UNIT=iw,FMT="(T2,A,T35,A,T80,A)") "-","T I M I N G","-"
       CASE DEFAULT
           STOP "BUG"
       END SELECT
       WRITE (UNIT=iw,FMT="(T2,A,T80,A)") "-","-"
       WRITE (UNIT=iw,FMT="(T2,A)") REPEAT("-",79)

       SELECT CASE(accu_type)
       CASE(accu_type_energy)
         WRITE (UNIT=iw,FMT="(T2,A,T35,A,T41,A,T45,2A18)")&
              "SUBROUTINE","CALLS"," ASD","SELF ENER","TOTAL ENER"
       CASE(accu_type_time)
         WRITE (UNIT=iw,FMT="(T2,A,T35,A,T41,A,T45,2A18)")&
              "SUBROUTINE","CALLS"," ASD","SELF TIME","TOTAL TIME"
       CASE DEFAULT
           STOP "BUG"
       END SELECT
       WRITE (UNIT=iw,FMT="(T45,4A9)")&
            "AVERAGE", "MAXIMUM", "AVERAGE","MAXIMUM"

       ! adjust fmt dynamically based on the max walltime.
       ! few clocks have more than 3 digits resolution, so stop there
       IF (num_routines > 0) THEN
         i = num_routines
         decimals = 3
         IF (max_walltime(i) >= 10000)   decimals = 2
         IF (max_walltime(i) >= 100000)  decimals = 1
         IF (max_walltime(i) >= 1000000) decimals = 0
         WRITE (UNIT=fmt,FMT="(A,I0,A)")&
           "(T2,A31,I7,1X,F4.1,4(1X,F8.",decimals,"))"
       END IF

       DO i=num_routines,1,-1

          j = indices(i)

          IF (max_walltime(i) >= mintime) THEN
             r_stat => list_get(timer_env%routine_stats, j)
             ! average stack depth
             asd = REAL(r_stat%stackdepth_accu, KIND=dp)/ &
                   REAL(MAX(1,r_stat%total_calls), KIND=dp)

             IF (sort_by_self_time) THEN
                WRITE (UNIT=iw,FMT=fmt) &
                 ADJUSTL(r_stat%routineN(1:31)),&
                 r_stat%total_calls,&
                 ASD,&
                 sum_ewalltime(j)/para_env%num_pe,&
                 max_ewalltime(i),&
                 sum_iwalltime(j)/para_env%num_pe,&
                 max_iwalltime(j)
             ELSE
                WRITE (UNIT=iw,FMT=fmt) &
                 ADJUSTL(r_stat%routineN(1:31)),&
                 r_stat%total_calls,&
                 ASD,&
                 sum_ewalltime(j)/para_env%num_pe,&
                 max_ewalltime(j),&
                 sum_iwalltime(j)/para_env%num_pe,&
                 max_iwalltime(i)
             END IF
          END IF

       END DO
       WRITE (UNIT=iw,FMT="(T2,A,/)") REPEAT("-",79)
       DEALLOCATE(indices)
    END IF

    NULLIFY (max_walltime)

  END SUBROUTINE timings_report_print

! *****************************************************************************
!> \brief Write accumulated callgraph information as cachegrind-file.
!> http://kcachegrind.sourceforge.net/cgi-bin/show.cgi/KcacheGrindCalltreeFormat
!> \par History
!>     12.2012  intial version[ole]
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE timings_report_callgraph(filename)
    CHARACTER(len=*), INTENT(in)             :: filename

    INTEGER, PARAMETER                       :: E = 1000, T = 100000

    INTEGER                                  :: i, unit
    TYPE(call_stat_type), POINTER            :: c_stat
    TYPE(dict_i4tuple_callstat_item_type), &
      DIMENSION(:), POINTER                  :: ct_items
    TYPE(routine_stat_type), POINTER         :: r_stat
    TYPE(timer_env_type), POINTER            :: timer_env

    CALL open_file(file_name=filename, file_status="REPLACE", file_action="WRITE", &
       file_form="FORMATTED", unit_number=unit)
    timer_env => get_timer_env()

    ! use outermost routine as total runtime
    r_stat => list_get(timer_env%routine_stats, 1)
    WRITE (UNIT=unit,FMT="(A)") "events: Walltime Energy"
    WRITE (UNIT=unit,FMT="(A,I0,1X,I0)") "summary: ", &
      INT(T*r_stat%incl_walltime_accu,KIND=int_8), &
      INT(E*r_stat%incl_energy_accu,KIND=int_8)

    DO i=1, list_size(timer_env%routine_stats)
       r_stat => list_get(timer_env%routine_stats, i)
       WRITE (UNIT=unit,FMT="(A,I0,A,A)") "fn=(",r_stat%routine_id,") ", r_stat%routineN
       WRITE (UNIT=unit,FMT="(A,I0,1X,I0)") "1 ", &
         INT(T*r_stat%excl_walltime_accu,KIND=int_8),&
         INT(E*r_stat%excl_energy_accu,KIND=int_8)
    END DO

    ct_items => dict_items(timer_env%callgraph)
    DO i=1, SIZE(ct_items)
       c_stat => ct_items(i)%value
       WRITE (UNIT=unit,FMT="(A,I0,A)") "fn=(", ct_items(i)%key(1), ")"
       WRITE (UNIT=unit,FMT="(A,I0,A)") "cfn=(", ct_items(i)%key(2), ")"
       WRITE (UNIT=unit,FMT="(A,I0,A)") "calls=", c_stat%total_calls, " 1"
       WRITE (UNIT=unit,FMT="(A,I0,1X,I0)") "1 ", &
         INT(T*c_stat%incl_walltime_accu,KIND=int_8),&
         INT(E*c_stat%incl_energy_accu,KIND=int_8)
    END DO
    DEALLOCATE(ct_items)

    CALL close_file(unit_number=unit, file_status="KEEP")

 END SUBROUTINE timings_report_callgraph
END MODULE timings_report

