!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Timing routines for accounting
!> \par History
!>      02.2004 made a stacked version (of stacks...) [Joost VandeVondele]
!>      11.2004 storable timer_envs (for f77 interface) [fawzi]
!>      10.2005 binary search to speed up lookup in timeset [fawzi]
!> \author JGH
! *****************************************************************************
MODULE timings
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_cputime,&
                                             m_memory,&
                                             m_walltime
  USE mathconstants,                   ONLY: zero
  USE message_passing,                 ONLY: mp_bcast
  USE string_utilities,                ONLY: uppercase
  USE util,                            ONLY: sort

  PRIVATE

  INTERFACE timeset
     MODULE PROCEDURE timeset_full, timeset_short
  END INTERFACE

  INTERFACE timestop
     MODULE PROCEDURE timestop_full, timestop_short
  END INTERFACE

  PUBLIC :: timeset, timestop, timeprint, trace_debug, print_stack,&
       use_HPM

  ! used to create a fresh timer_env
  PUBLIC :: timer_env_type, timer_env_p_type
  PUBLIC :: add_timer_env, rm_timer_env, get_timer_env
  PUBLIC :: timer_env_create, timer_env_retain, timer_env_release

  ! want to use HPM if available ? needs care switching true/false (i.e. use_HPM=.true. timeset ... timestop use_HPM=false)
  LOGICAL :: use_HPM = .FALSE.

  ! hard coded stack sizes
  INTEGER, PARAMETER :: max_timer_env_stack=10  ! number of stacked timers
  INTEGER, PARAMETER :: max_stack = 50          ! the calling stack max depth for a given timer env
  INTEGER, PARAMETER :: max_timer = 2000        ! total number of different callers

! *****************************************************************************
  TYPE timer_env_type
     !PRIVATE
     INTEGER :: ref_count, id_nr
     ! subroutine calling stack
     INTEGER :: stack_size
     CHARACTER ( LEN = 80 ) :: routine_stack ( max_stack )
     INTEGER :: handle_stack (max_stack)
     ! subroutine timing and performance information
     CHARACTER ( LEN = 80 ) :: clock_name ( max_timer )
     INTEGER :: sort_index ( max_timer )
     INTEGER :: routine_calls ( max_timer )
     INTEGER :: routine_stack_depth ( max_timer )
     INTEGER :: clock_status ( max_timer )
     REAL (KIND=dp) :: clock_start ( max_timer )
     REAL (KIND=dp) :: clock_accu ( max_timer )
     REAL (KIND=dp) :: wallclock_start ( max_timer )
     REAL (KIND=dp) :: wallclock_accu ( max_timer )
     INTEGER :: clock_type ( max_timer )
     REAL (KIND=dp) :: perf_count ( max_timer )
     CHARACTER ( LEN = 10 ) :: perf_unit ( max_timer )
     INTEGER :: num_timer 
     CHARACTER ( LEN = 15 ) :: astring, bstring
     LOGICAL :: trace,trace_sync
     INTEGER :: mepos,group,num_pe
  END TYPE timer_env_type

! *****************************************************************************
  TYPE timer_env_p_type
     TYPE(timer_env_type), POINTER :: timer_env
  END TYPE timer_env_p_type

  INTEGER, SAVE, PRIVATE :: stack_pointer = 0       ! position in timer_env_stack
  INTEGER, SAVE, PRIVATE :: last_timer_env_id = 0
  TYPE(timer_env_p_type), DIMENSION(max_timer_env_stack) :: timer_env_stack 
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'timings'

CONTAINS

! *****************************************************************************
!> \brief adds the given timer_env to the top of the stack
!> \note
!>      for each init_timer_env there should be the symmetric call to 
!>      rm_timer_env
!> \par History
!>      02.2004 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE add_timer_env(timer_env)
    TYPE(timer_env_type), OPTIONAL, POINTER  :: timer_env

    stack_pointer=stack_pointer+1
    IF (stack_pointer>max_timer_env_stack) THEN
       PRINT*,'stack_pointer too large is add_timer_env'
       STOP 1
    ENDIF
    NULLIFY(timer_env_stack(stack_pointer)%timer_env)
    IF (PRESENT(timer_env)) THEN
       timer_env_stack(stack_pointer)%timer_env => timer_env
       IF (ASSOCIATED(timer_env)) THEN
          CALL timer_env_retain(timer_env)
       END IF
    END IF
    IF (.NOT.ASSOCIATED(timer_env_stack(stack_pointer)%timer_env)) THEN
       CALL timer_env_create(timer_env_stack(stack_pointer)%timer_env)
    END IF
  END SUBROUTINE add_timer_env

! *****************************************************************************
!> \brief creates a new timer env
!> \author fawzi
! *****************************************************************************
  SUBROUTINE timer_env_create(timer_env)
    TYPE(timer_env_type), POINTER            :: timer_env

    INTEGER                                  :: stat

    ALLOCATE(timer_env,stat=stat)
    IF (stat/=0) THEN
       PRINT*,moduleN//':timer_env_create error allocationg timer_env'
       STOP 1
    END IF
    last_timer_env_id=last_timer_env_id+1
    timer_env%id_nr=last_timer_env_id
    timer_env%ref_count=1
    timer_env%stack_size=0
    timer_env%handle_stack=0
    timer_env%num_timer=0
    timer_env%trace=.FALSE. 
    timer_env%trace_sync=.FALSE. 
  END SUBROUTINE timer_env_create

! *****************************************************************************
!> \brief removes the current timer env from the stack
!> \note
!>      for each rm_timer_env there should have been the symmetric call to 
!>      add_timer_env
!> \par History
!>      02.2004 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE rm_timer_env()
    IF (stack_pointer<1) THEN
       PRINT*,'no timer env available in destroy_timer_env'
       STOP 1
    ENDIF
    CALL timer_env_release(timer_env_stack(stack_pointer)%timer_env)
    stack_pointer=stack_pointer-1
  END SUBROUTINE rm_timer_env

! *****************************************************************************
!> \brief returns the current timer env from the stack
!> \author fawzi
! *****************************************************************************
  FUNCTION get_timer_env() RESULT(res)
    TYPE(timer_env_type), POINTER            :: res

    IF (stack_pointer<1) THEN
       PRINT*,'no timer env available in get_timer_env'
       STOP 1
    ENDIF
    res => timer_env_stack(stack_pointer)%timer_env
  END FUNCTION get_timer_env

! *****************************************************************************
!> \brief retains the given timer env
!> \param timer_env the timer env to retain
!> \author fawzi
! *****************************************************************************
  SUBROUTINE timer_env_retain(timer_env)
    TYPE(timer_env_type), POINTER            :: timer_env

    IF (.NOT.ASSOCIATED(timer_env)) THEN
       PRINT*,moduleN,':timer_env_retain ERROR, unassociated timer_env'
       STOP 1
    ENDIF
    IF (timer_env%ref_count<1) THEN
       PRINT*,moduleN,':timer_env_retain ERROR, timer_env%ref_count=',&
            timer_env%ref_count
       STOP 1
    ENDIF
    timer_env%ref_count=timer_env%ref_count+1
  END SUBROUTINE timer_env_retain

! *****************************************************************************
!> \brief releases the given timer env
!> \param timer_env the timer env to release
!> \author fawzi
! *****************************************************************************
  SUBROUTINE timer_env_release(timer_env)
    TYPE(timer_env_type), POINTER            :: timer_env

    INTEGER                                  :: stat

    IF (ASSOCIATED(timer_env)) THEN
       IF (timer_env%ref_count<1) THEN
          PRINT*,moduleN,':timer_env_release ERROR, timer_env%ref_count=',&
               timer_env%ref_count,"id_nr=",timer_env%id_nr
          STOP 1
       ENDIF
       timer_env%ref_count=timer_env%ref_count-1
       IF (timer_env%ref_count==0) THEN
          timer_env%stack_size=0
          timer_env%handle_stack=0
          timer_env%num_timer=0
          timer_env%trace=.FALSE.
          timer_env%trace_sync=.FALSE.
          DEALLOCATE(timer_env,stat=stat)
          IF (stat/=0) THEN
             PRINT*,moduleN,':timer_env_release ERROR deallocating, stat=',stat
             STOP 1
          ENDIF
       END IF
    END IF
    NULLIFY(timer_env)
  END SUBROUTINE timer_env_release

! *****************************************************************************
!> \brief Start timer (short call). Provides the usual shortcut
!> \par History
!>      none
!> \author JGH
! *****************************************************************************
  SUBROUTINE timeset_short(name, handle )
    CHARACTER(LEN=*), INTENT(IN)             :: name
    INTEGER                                  :: handle

    CALL timeset_full(name, "I", "", handle)
  END SUBROUTINE timeset_short

! *****************************************************************************
!> \brief Start timer (long call)
!> \par History
!>      none
!> \author JGH
! *****************************************************************************
  SUBROUTINE timeset_full ( name, ct, pu, handle )
    CHARACTER(LEN=*), INTENT(IN)             :: name
    CHARACTER(LEN=1), INTENT(IN)             :: ct           ! Clock type [0|1]
    CHARACTER(LEN=*), INTENT(IN)             :: pu           ! Performance unit
    INTEGER, INTENT(OUT)                     :: handle

    CHARACTER(LEN=20)                        :: mytag
    CHARACTER(LEN=40)                        :: sformat
    CHARACTER(LEN=400)                       :: line, line2, mystring
    INTEGER                                  :: i, stack_size
    TYPE(timer_env_type), POINTER            :: timer_env
#if defined(__HPM)
#include "f_hpm.h"
#endif

    IF (stack_pointer<1) THEN
       PRINT*,'no timer env available in timeset'
       STOP 1
    ENDIF
    IF (LEN_TRIM(name)>LEN(timer_env%clock_name)) THEN
       PRINT *,"timings:timeset, timer name too long: '",TRIM(name),"'"
       STOP 1
    END IF
    timer_env => timer_env_stack(stack_pointer)%timer_env
    !..add name to subroutine stack
    timer_env%stack_size = timer_env%stack_size + 1
    stack_size = timer_env%stack_size
    timer_env%routine_stack ( stack_size ) = name

    !..binary search of the timer
    handle = timer_binsearch(name,timer_env%clock_name,timer_env%sort_index,&
         1,timer_env%num_timer)
    IF (handle>0) THEN
       IF (timer_env%clock_name(timer_env%sort_index(handle))==name) go to 33
    END IF
    !..it was a new timer
    IF (timer_env%num_timer>=max_timer) THEN
       PRINT *,'ERROR hit max_timer in timeset'
       PRINT *,"name ",name,timer_env%num_timer
       PRINT *,timer_env%clock_name
       STOP 1
    ENDIF

    DO i=timer_env%num_timer,handle+1,-1
       timer_env%sort_index(i+1)=timer_env%sort_index(i)
    END DO
    timer_env%sort_index(handle+1)=timer_env%num_timer+1
    handle=handle+1

    timer_env%num_timer = timer_env%num_timer + 1
    timer_env%clock_name ( timer_env%num_timer ) = &
         name
    timer_env%perf_unit ( timer_env%num_timer ) = pu
    IF ( ct == 'I' ) timer_env%clock_type ( timer_env%num_timer ) = 0
    IF ( ct == 'E' ) timer_env%clock_type ( timer_env%num_timer ) = 1
    timer_env%clock_status ( timer_env%num_timer ) = 0
    timer_env%clock_accu ( timer_env%num_timer ) = zero
    timer_env%wallclock_accu ( timer_env%num_timer ) = zero
    timer_env%perf_count ( timer_env%num_timer ) = zero
    timer_env%routine_calls ( timer_env%num_timer ) = 0
    timer_env%routine_stack_depth ( timer_env%num_timer ) = 0

33  handle=timer_env%sort_index(handle)

#if defined(__HPM)
    IF (use_HPM) THEN
       mytag=TRIM(name)
       mytag(20:20)=CHAR(0)
       CALL f_hpmstart(handle,mytag)
    ENDIF
#endif

    !..start clock (if not yet started, i.e. recursive calls)
    IF ( timer_env%clock_status ( handle ) ==0 ) THEN
       timer_env%clock_start ( handle ) = m_cputime ( )
       timer_env%wallclock_start ( handle ) = m_walltime ( )
       timer_env%clock_status ( handle ) = 1
    END IF

    timer_env%routine_calls ( handle ) = &
         timer_env%routine_calls ( handle ) + 1
    timer_env%routine_stack_depth ( handle ) =&
         timer_env%routine_stack_depth ( handle ) + stack_size

    timer_env%handle_stack (stack_size) = handle

    !..if debug mode echo the subroutine name
    IF ( timer_env%trace ) THEN
       WRITE(sformat,*) "(A,",MAX(1,3*stack_size-9),"X,I4,1X,I6,1X,A,A)"
       WRITE ( mystring, sformat) &
            timer_env%astring,stack_size, &
            timer_env%routine_calls ( handle ), &
            TRIM( timer_env%routine_stack ( stack_size )),"       start"
       WRITE(line,'(A,1X,I0,1X,A2)') TRIM(mystring),(m_memory()+1023)/1024,"Kb"
       IF (timer_env%trace_sync) THEN
          DO i=0,timer_env%num_pe-1
             line2=line
             CALL mp_bcast(line2,i,timer_env%group)
             IF (timer_env%mepos==0) THEN
                WRITE(6,*) TRIM(line2)
             ENDIF
          ENDDO
       ELSE
          WRITE(6,*) TRIM(line)
       ENDIF
    ENDIF

  END SUBROUTINE timeset_full

! *****************************************************************************
!> \brief End timer
!> \par History
!>      none
!> \author JGH
! *****************************************************************************
  SUBROUTINE timestop_short( handle ) ! provide the usual shortcut
    INTEGER                                  :: handle

    CALL timestop_full(0.0_dp,handle)
  END SUBROUTINE timestop_short

! *****************************************************************************
!> \brief End timer (long call)
!> \par History
!>      none
!> \author JGH
! *****************************************************************************
  SUBROUTINE timestop_full ( perf, handle )
    REAL(KIND=dp)                            :: perf
    INTEGER                                  :: handle

    CHARACTER(LEN=40)                        :: sformat
    CHARACTER(LEN=400)                       :: line, line2, mystring
    INTEGER                                  :: i
    REAL(KIND=dp)                            :: cpt, wt
    TYPE(timer_env_type), POINTER            :: timer_env
#if defined(__HPM)
#include "f_hpm.h"
#endif

    IF (stack_pointer<1) THEN
       PRINT*,'no timer env available in timestop'
       STOP  1
    ENDIF

    timer_env => timer_env_stack(stack_pointer)%timer_env

    IF ( handle /= timer_env%handle_stack (timer_env%stack_size) ) THEN
       PRINT *, "timer_env%stack_size",timer_env%stack_size," handle ",handle," timer_env%id_nr ",timer_env%id_nr
       PRINT*,'mismatched timestop '&
            //TRIM(timer_env%routine_stack(timer_env%stack_size))//&
            ' in routine timestop'
       STOP 1
    END IF
    timer_env%handle_stack (timer_env%stack_size)=0

    IF ( handle > timer_env%num_timer ) THEN
       ! intermediate work-around
       PRINT*,'handle > num_timer in routine timestop'
       STOP 1
    ELSE
#if defined(__HPM)
       IF (use_HPM) THEN
          CALL f_hpmstop(handle)
       ENDIF
#endif
       cpt = m_cputime ( ) - timer_env%clock_start ( handle )
       wt = m_walltime ( ) - timer_env%wallclock_start ( handle )
       timer_env%clock_accu ( handle ) = &
            timer_env%clock_accu ( handle ) + cpt
       timer_env%wallclock_accu ( handle ) = &
            timer_env%wallclock_accu ( handle ) + wt
       timer_env%clock_status ( handle ) = 0
       timer_env%perf_count ( handle ) = &
            timer_env%perf_count ( handle ) + perf

       DO i = 1, timer_env% num_timer
          IF ( timer_env%clock_status ( i ) == 1 .AND. &
               timer_env%clock_type ( i ) == 1 ) THEN
             timer_env%clock_start ( i ) = &
                  timer_env%clock_start ( i ) + cpt
             timer_env%wallclock_start ( i ) = &
                  timer_env%wallclock_start ( i ) + wt
          END IF
       END DO

       !..if debug mode echo the subroutine name
       IF ( timer_env%trace ) THEN
          WRITE(sformat,*) "(A,",MAX(1,3*timer_env%stack_size-9),"X,I4,1X,I6,1X,A,F12.3)"
          WRITE ( mystring, sformat) &
               timer_env%bstring, & 
               timer_env%stack_size, &
               timer_env%routine_calls ( handle ), &
               TRIM(timer_env%routine_stack(timer_env%stack_size)),cpt
          WRITE(line,'(A,1X,I0,1X,A2)') TRIM(mystring),(m_memory()+1023)/1024,"Kb"
          IF (timer_env%trace_sync) THEN
             DO i=0,timer_env%num_pe-1
                line2=line
                CALL mp_bcast(line2,i,timer_env%group)
                IF (timer_env%mepos==0) THEN
                   WRITE(6,*) TRIM(line2)
                ENDIF
             ENDDO
          ELSE
             WRITE(6,*) TRIM(line)
          ENDIF
       ENDIF

       !..remove name from subroutine stack
       timer_env%stack_size = timer_env%stack_size - 1
    END IF
  END SUBROUTINE timestop_full

! *****************************************************************************
!> \brief Print accumulated information on timers
!> \par History
!>      none
!> \author JGH
! *****************************************************************************
  SUBROUTINE timeprint(iw,r_timings)
    INTEGER, INTENT(IN)                      :: iw
    REAL(KIND=dp), INTENT(IN)                :: r_timings

    INTEGER                                  :: i, j
    CHARACTER(LEN=1), DIMENSION(0:1)         :: ct = (/"I","E"/)
    INTEGER, DIMENSION(max_timer)            :: index
    REAL(KIND=dp)                            :: asd, maxtime, mintime, perf
    TYPE(timer_env_type), POINTER            :: timer_env

   IF (stack_pointer<1) THEN
       PRINT*,'no timer env available in time_set'
       STOP 1
    ENDIF

    timer_env => timer_env_stack(stack_pointer)%timer_env

    IF (timer_env%num_timer > 0) THEN

       CALL sort(timer_env%wallclock_accu,timer_env%num_timer,index)

       maxtime = timer_env%wallclock_accu(timer_env%num_timer)
       mintime = maxtime*r_timings

       WRITE (UNIT=iw,FMT="(/,T2,A)") REPEAT("-",79)
       WRITE (UNIT=iw,FMT="(T2,A,T80,A)") "-","-"
       WRITE (UNIT=iw,FMT="(T2,A,T35,A,T80,A)") "-","T I M I N G","-"
       WRITE (UNIT=iw,FMT="(T2,A,T80,A)") "-","-"
       WRITE (UNIT=iw,FMT="(T2,A)") REPEAT("-",79)
       WRITE (UNIT=iw,FMT="(T2,A,T32,A,T38,A,T45,A,T57,A,T69,A)")&
            "SUBROUTINE","CALLS"," ASD","CPU TIME T","PERFORMANCE","ELAPSED TIME"

       DO i=timer_env%num_timer,1,-1

          j = INDEX(i)

          IF (timer_env%wallclock_accu(i) > mintime) THEN

             IF (timer_env%wallclock_accu(i) > zero) THEN
                perf = timer_env%perf_count(j)/timer_env%wallclock_accu(i)
             ELSE
                perf = zero
             END IF
             ! average stack depth
             asd=REAL(timer_env%routine_stack_depth(j),KIND=dp)/ &
                  REAL(MAX(1,timer_env%routine_calls(j)),KIND=dp)

             IF (LEN_TRIM(timer_env%perf_unit(j)) == 0) THEN
                WRITE (UNIT=iw,&
                     FMT="(T2,A,T30,I7,T38,F4.1,T42,F11.2,T54,A,T62,A,T70,F11.2)")&
                     TRIM(timer_env%clock_name(j)),timer_env%routine_calls(j), & 
                     ASD,timer_env%clock_accu(j),&
                     ct(timer_env%clock_type(j)),"---",timer_env%wallclock_accu(i)
             ELSE
                WRITE (UNIT=iw,&
                     FMT="(T2,A,T30,I7,T38,F4.1,T42,F11.2,T54,A,T56,F5.1,T62,A,T70,F11.2)")&
                     TRIM(timer_env%clock_name(j)),timer_env%routine_calls(j),&
                     ASD,timer_env%clock_accu(j),&
                     ct(timer_env%clock_type(j)),perf,TRIM(timer_env%perf_unit(j)),&
                     timer_env%wallclock_accu(i)
             END IF
          END IF

       END DO
       WRITE (UNIT=iw,FMT="(T2,A,/)") REPEAT("-",79)
    END IF

  END SUBROUTINE timeprint

! *****************************************************************************
!> \brief Set routine tracer
!> \par History
!>      none
!> \author JGH
! *****************************************************************************
  SUBROUTINE trace_debug ( set, trace_sync, para_env )
    CHARACTER(len=*)                         :: set
    LOGICAL                                  :: trace_sync
    TYPE(cp_para_env_type)                   :: para_env

    CHARACTER(len=5)                         :: set2
    TYPE(timer_env_type), POINTER            :: timer_env

    IF (stack_pointer<1) THEN
       PRINT*,'no timer env available in trace_debug'
       STOP 1
    ENDIF

    timer_env => timer_env_stack(stack_pointer)%timer_env

    set2=TRIM(set)
    CALL uppercase(set2)

    timer_env%trace_sync=trace_sync

    SELECT CASE (set2)
    CASE ("START")
       timer_env%trace = .TRUE.
    CASE ("STOP")
       timer_env%trace = .FALSE.
    CASE DEFAULT
       ! intermediate work-around
       PRINT*,'unknown action requested in routine trace_debug'
       STOP 1
    END SELECT

    timer_env%astring = " >>>>>>>>>>>>> "
    timer_env%bstring = " <<<<<<<<<<<<< "
    WRITE (timer_env%astring(1:9),"(I4,A,I4)") para_env%mepos,":",&
         para_env%num_pe
    WRITE (timer_env%bstring(1:9),"(I4,A,I4)") para_env%mepos,":",&
         para_env%num_pe

    timer_env%mepos=para_env%mepos
    timer_env%group=para_env%group
    timer_env%num_pe=para_env%num_pe

  END SUBROUTINE trace_debug

! *****************************************************************************
!> \brief Print current routine stack
!> \par History
!>      none
!> \author JGH
! *****************************************************************************
  SUBROUTINE print_stack ( iw )

    INTEGER                                  :: iw

    INTEGER                                  :: i, io
    TYPE(timer_env_type), POINTER            :: timer_env

    io = iw
    IF (stack_pointer<1) THEN
       PRINT*,'no timer env available in print_stack'
       STOP 1
    ENDIF

    timer_env => timer_env_stack(stack_pointer)%timer_env

    WRITE ( io , '(/,A,/)') " ===== Routine Calling Stack ===== "
    DO i = timer_env%stack_size, 1, -1
       WRITE ( io, '(T10,I4,1X,A)') i,TRIM(timer_env%routine_stack(i))
    END DO

  END SUBROUTINE print_stack

! *****************************************************************************
!> \brief makes a binary search in an array of strings ordered with the index
!>      diven in ord_idx and returns the index of arr_idx
!>      where the element is (or the element after which it should be inserted)
!> \param el the element to search
!> \param list the list of strings
!> \param ord_idx indexes that make list ordered in growing order
!> \param lb lower bound for the search (defaults to 1)
!> \param ub upper bound for the search (defaults to size(list))
!> \author fawzi
! *****************************************************************************
  FUNCTION timer_binsearch(el,list,ord_idx,lb,ub)RESULT(res)
    CHARACTER(len=*), INTENT(in)             :: el
    CHARACTER(len=80), DIMENSION(:), &
      INTENT(in)                             :: list
    INTEGER, DIMENSION(:), INTENT(in)        :: ord_idx
    INTEGER, INTENT(in)                      :: lb, ub
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'timer_binsearch', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: li, mi, ui
    LOGICAL                                  :: failure

    failure=.FALSE.
    li=lb
    ui=ub
    IF (ui>SIZE(list)) THEN
       PRINT *, "error in timings:timer_binsearch, ui>size(list)",ui,SIZE(list)
       STOP 1
    END IF
    DO
       IF (li>ui) EXIT
       mi=(li+ui)/2
       IF (el<list(ord_idx(mi))) THEN
          ui=mi-1
       ELSE
          li=mi+1
       END IF
    END DO
    res=ui
  END FUNCTION timer_binsearch

END MODULE timings

