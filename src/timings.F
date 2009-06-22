!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
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
  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE machine,                         ONLY: m_flush,&
                                             m_memory,&
                                             m_walltime
  USE mathconstants,                   ONLY: zero
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum
  USE string_utilities,                ONLY: uppercase
  USE util,                            ONLY: sort

  PRIVATE

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
  INTEGER, PARAMETER :: max_clock_name = 80     ! length of strings, needs adjusting of formats as well

! *****************************************************************************
  TYPE timer_env_type
     PRIVATE
     INTEGER                                         :: ref_count, id_nr
     ! subroutine calling stack
     INTEGER                                         :: stack_size
     CHARACTER ( LEN = max_clock_name )              :: routine_stack (max_stack)
     REAL (KIND=dp)                                  :: wallclock_start (max_stack)
     INTEGER                                         :: handle_stack (max_stack)
     ! timing info
     CHARACTER ( LEN = max_clock_name )              :: clock_name (max_timer)
     INTEGER                                         :: sort_index ( max_timer )
     INTEGER                                         :: routine_calls ( max_timer )
     INTEGER                                         :: routine_stack_depth ( max_timer )
     INTEGER                                         :: clock_status ( max_timer )
     REAL (KIND=dp)                                  :: exclusive_wallclock_accu ( max_timer )
     REAL (KIND=dp)                                  :: inclusive_wallclock_accu ( max_timer )
     INTEGER                                         :: num_timer 
     CHARACTER ( LEN = 15 )                          :: astring, bstring
     LOGICAL                                         :: trace,trace_sync,trace_master
     INTEGER                                         :: mepos,group,num_pe
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
    timer_env%trace_master=.FALSE. 
    timer_env%group=-1
    timer_env%mepos=-1
    timer_env%num_pe=-1
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
          timer_env%trace_master=.FALSE.
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
!> \brief Start timer
!> \par History
!>      none
!> \author JGH
! *****************************************************************************
  SUBROUTINE timeset ( name, handle )
    CHARACTER(LEN=*), INTENT(IN)             :: name
    INTEGER, INTENT(OUT)                     :: handle

    CHARACTER(LEN=20)                        :: mytag
    CHARACTER(LEN=40)                        :: sformat
    CHARACTER(LEN=400)                       :: line, line2, mystring
    INTEGER                                  :: i, stack_size
    INTEGER(KIND=int_8)                      :: mem
    TYPE(timer_env_type), POINTER            :: timer_env
    REAL(KIND=dp)                            :: wt_now
#if defined(__HPM)
#include "f_hpm.h"
#endif

    wt_now=m_walltime ()

    IF (stack_pointer<1) THEN
       PRINT*,'no timer env available in timeset'
       STOP 1
    ENDIF
    IF (LEN_TRIM(name)>LEN(timer_env%clock_name)) THEN
       PRINT *,"timings:timeset, timer name too long: '",TRIM(name),"'"
       STOP 1
    END IF
    timer_env => timer_env_stack(stack_pointer)%timer_env

    !..binary search of the timer
    handle = timer_binsearch(name,timer_env%clock_name,timer_env%sort_index,&
         1,timer_env%num_timer)
    IF (handle>0) THEN
       IF (timer_env%clock_name(timer_env%sort_index(handle))==name) go to 33
    END IF

    !.. initialize the newly found clock_name
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
    timer_env%clock_name ( timer_env%num_timer ) = name
    timer_env%clock_status ( timer_env%num_timer ) = 0
    timer_env%exclusive_wallclock_accu ( timer_env%num_timer ) = zero
    timer_env%inclusive_wallclock_accu ( timer_env%num_timer ) = zero
    timer_env%routine_calls ( timer_env%num_timer ) = 0
    timer_env%routine_stack_depth ( timer_env%num_timer ) = 0

    ! this is the handle of the current clock
33  handle=timer_env%sort_index(handle)

#if defined(__HPM)
    IF (use_HPM) THEN
       mytag=TRIM(name)
       mytag(20:20)=CHAR(0)
       CALL f_hpmstart(handle,mytag)
    ENDIF
#endif

    !..add name to subroutine stack
    timer_env%stack_size = timer_env%stack_size + 1
    stack_size = timer_env%stack_size
    timer_env%routine_stack ( stack_size ) = name
    timer_env%wallclock_start ( stack_size ) = wt_now
    timer_env%handle_stack (stack_size) = handle

    ! count the number of active instances (recursion)
    timer_env%clock_status ( handle ) = timer_env%clock_status ( handle ) + 1
   
    ! some properties
    timer_env%routine_calls ( handle ) = timer_env%routine_calls ( handle ) + 1
    timer_env%routine_stack_depth ( handle ) = timer_env%routine_stack_depth ( handle ) + stack_size

    IF (timer_env%routine_stack(stack_size).NE. timer_env%clock_name(timer_env%handle_stack (stack_size))) THEN
       STOP "CP2K BUG with timeset/timestop?"
    ENDIF

    !..if debug mode echo the subroutine name
    IF ( timer_env%trace ) THEN
       WRITE(sformat,*) "(A,",MAX(1,3*stack_size-9),"X,I4,1X,I6,1X,A,A)"
       WRITE ( mystring, sformat) &
            timer_env%astring,stack_size, &
            timer_env%routine_calls ( handle ), &
            TRIM( timer_env%routine_stack ( stack_size )),"       start"
       mem=m_memory()
       WRITE(line,'(A,1X,I0,1X,A2)') TRIM(mystring),(mem+1023)/1024,"Kb"
       IF (timer_env%trace_sync) THEN
          DO i=0,timer_env%num_pe-1
             line2=line
             CALL mp_bcast(line2,i,timer_env%group)
             IF (timer_env%mepos==0) THEN
                WRITE(6,*) TRIM(line2)
                CALL m_flush(6)
             ENDIF
          ENDDO
       ELSE
          IF (.NOT.timer_env%trace_master .OR. timer_env%mepos==0) THEN
             WRITE(6,*) TRIM(line)
             CALL m_flush(6)
          ENDIF
       ENDIF
    ENDIF

  END SUBROUTINE timeset

! *****************************************************************************
!> \brief End timer
!> \par History
!>      none
!> \author JGH
! *****************************************************************************
  SUBROUTINE timestop( handle )
    INTEGER                                  :: handle

    CHARACTER(LEN=40)                        :: sformat
    CHARACTER(LEN=400)                       :: line, line2, mystring
    INTEGER                                  :: handle_up,i
    INTEGER(KIND=int_8)                      :: mem
    REAL(KIND=dp)                            :: elapsed_wt, wt_now
    TYPE(timer_env_type), POINTER            :: timer_env

#if defined(__HPM)
#include "f_hpm.h"
#endif

    wt_now = m_walltime()

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

    IF ( handle > timer_env%num_timer ) THEN
       ! intermediate work-around
       PRINT*,'handle > num_timer in routine timestop'
       STOP 1
    ENDIF

#if defined(__HPM)
    IF (use_HPM) THEN
       CALL f_hpmstop(handle)
    ENDIF
#endif

    ! add the elapsed time for this timeset/timestop to the time accumulator
    elapsed_wt = wt_now - timer_env%wallclock_start ( timer_env%stack_size )

    ! if we're the last instance in the stack, we do the accounting of the total time
    IF (timer_env%clock_status ( handle ) == 1) THEN
       timer_env%inclusive_wallclock_accu ( handle ) = timer_env%inclusive_wallclock_accu ( handle ) + elapsed_wt
    ENDIF

    ! exclusive time we always sum, since children will correct this time with their total time
    timer_env%exclusive_wallclock_accu ( handle ) = timer_env%exclusive_wallclock_accu ( handle ) + elapsed_wt

    ! we fixup the clock of the caller
    IF (timer_env%stack_size>1) THEN
       handle_up=timer_env%handle_stack (timer_env%stack_size-1) 
       timer_env%exclusive_wallclock_accu ( handle_up ) = timer_env%exclusive_wallclock_accu ( handle_up ) - elapsed_wt
    ENDIF

    !..if debug mode echo the subroutine name
    IF ( timer_env%trace ) THEN
       WRITE(sformat,*) "(A,",MAX(1,3*timer_env%stack_size-9),"X,I4,1X,I6,1X,A,F12.3)"
       WRITE ( mystring, sformat) &
            timer_env%bstring, & 
            timer_env%stack_size, &
            timer_env%routine_calls ( handle ), &
            TRIM(timer_env%routine_stack(timer_env%stack_size)),elapsed_wt
       mem=m_memory()
       WRITE(line,'(A,1X,I0,1X,A2)') TRIM(mystring),(mem+1023)/1024,"Kb"
       IF (timer_env%trace_sync) THEN
          DO i=0,timer_env%num_pe-1
             line2=line
             CALL mp_bcast(line2,i,timer_env%group)
             IF (timer_env%mepos==0) THEN
                WRITE(6,*) TRIM(line2)
                CALL m_flush(6)
             ENDIF
          ENDDO
       ELSE
          IF (.NOT.timer_env%trace_master .OR. timer_env%mepos==0) THEN
             WRITE(6,*) TRIM(line)
             CALL m_flush(6)
          ENDIF
       ENDIF
    ENDIF

    !..remove name from subroutine stack
    timer_env%clock_status ( handle ) = timer_env%clock_status ( handle ) - 1
    timer_env%handle_stack (timer_env%stack_size)=0
    timer_env%stack_size = timer_env%stack_size - 1

  END SUBROUTINE timestop

! *****************************************************************************
!> \brief Print accumulated information on timers
!> \param para_env is needed here because the group in the timer_env is not valid
!> \par History
!>      none
!> \author JGH
! *****************************************************************************
  SUBROUTINE timeprint(iw,r_timings,para_env)
    INTEGER, INTENT(IN)                      :: iw
    REAL(KIND=dp), INTENT(IN)                :: r_timings
    TYPE(cp_para_env_type), INTENT(IN)       :: para_env

    CHARACTER(LEN=max_clock_name)            :: iclock
    INTEGER                                  :: i, j, num_timer
    INTEGER, DIMENSION(max_timer)            :: index
    REAL(KIND=dp)                            :: asd, imax_ewall, imax_iwall, &
                                                isum_ewall, isum_iwall, &
                                                maxtime, mintime
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: max_ewalltime, max_iwalltime, &
                                                sum_ewalltime, sum_iwalltime
    TYPE(timer_env_type), POINTER            :: timer_env

    IF (stack_pointer<1) THEN
       PRINT*,'no timer env available in time_set'
       STOP 1
    ENDIF

    timer_env => timer_env_stack(stack_pointer)%timer_env

    ! find the number of timing events on the parent
    ! note: we assume that timing events on the parent
    !       include all possible timing events.. i.e. there
    !       are no other timing events on the rest of the group.
    !       If this is the case they will be ignored in the 
    !       timing report.
    num_timer=timer_env%num_timer
    CALL mp_bcast(num_timer,para_env%source,para_env%group)
    
    ALLOCATE(max_iwalltime(num_timer))
    ALLOCATE(sum_iwalltime(num_timer))
    ALLOCATE(max_ewalltime(num_timer))
    ALLOCATE(sum_ewalltime(num_timer))

    ! for each clock_name find the needed quantities over all CPUs.
    DO i=1,num_timer
       iclock=timer_env%clock_name(i)
       CALL mp_bcast(iclock,para_env%source,para_env%group)
       imax_iwall=0.0_dp
       isum_iwall=0.0_dp
       imax_ewall=0.0_dp
       isum_ewall=0.0_dp
       DO j=1,timer_env%num_timer
          IF (timer_env%clock_name(j)==iclock) THEN
             imax_iwall=timer_env%inclusive_wallclock_accu(j)
             isum_iwall=timer_env%inclusive_wallclock_accu(j)
             imax_ewall=timer_env%exclusive_wallclock_accu(j)
             isum_ewall=timer_env%exclusive_wallclock_accu(j)
             EXIT
          ENDIF
       ENDDO
       CALL mp_sum(isum_iwall,para_env%group)
       CALL mp_max(imax_iwall,para_env%group)
       max_iwalltime(i)=imax_iwall
       sum_iwalltime(i)=isum_iwall

       CALL mp_sum(isum_ewall,para_env%group)
       CALL mp_max(imax_ewall,para_env%group)
       max_ewalltime(i)=imax_ewall
       sum_ewalltime(i)=isum_ewall
    ENDDO

    IF (num_timer > 0 .AND. iw>0) THEN

       CALL sort(max_iwalltime,num_timer,index)

       maxtime = max_iwalltime(num_timer)
       mintime = maxtime*r_timings

       WRITE (UNIT=iw,FMT="(/,T2,A)") REPEAT("-",79)
       WRITE (UNIT=iw,FMT="(T2,A,T80,A)") "-","-"
       WRITE (UNIT=iw,FMT="(T2,A,T35,A,T80,A)") "-","T I M I N G","-"
       WRITE (UNIT=iw,FMT="(T2,A,T80,A)") "-","-"
       WRITE (UNIT=iw,FMT="(T2,A)") REPEAT("-",79)

       WRITE (UNIT=iw,FMT="(T2,A,T34,A,T40,A,T45,2A18)")&
            "SUBROUTINE","CALLS"," ASD","SELF TIME","TOTAL TIME"
       WRITE (UNIT=iw,FMT="(T45,4A9)")&
            "AVERAGE", "MAXIMUM", "AVERAGE","MAXIMUM"

       DO i=num_timer,1,-1

          j = INDEX(i)

          IF (max_iwalltime(i) > mintime) THEN

             ! average stack depth
             asd=REAL(timer_env%routine_stack_depth(j),KIND=dp)/ &
                  REAL(MAX(1,timer_env%routine_calls(j)),KIND=dp)

             WRITE (UNIT=iw,FMT="(T2,A,T33,I6,T40,F4.1,T45,4F9.2)") &
                          TRIM(timer_env%clock_name(j)),&
                          timer_env%routine_calls(j),&
                          ASD,&
                          sum_ewalltime(j)/para_env%num_pe,&
                          max_ewalltime(j),&
                          sum_iwalltime(j)/para_env%num_pe,&
                          max_iwalltime(i)
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
  SUBROUTINE trace_debug ( set, trace_sync, trace_master, para_env )
    CHARACTER(len=*)                         :: set
    LOGICAL                                  :: trace_sync, trace_master
    TYPE(cp_para_env_type), INTENT(IN)       :: para_env

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
    timer_env%trace_master=trace_master
    ! this is a hack, timer_env can not retain para_env in a clean way
    ! as such, there is little guarantee this communicator remains valid
    timer_env%group=para_env%group
    timer_env%mepos=para_env%mepos
    timer_env%num_pe=para_env%num_pe

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
    WRITE (timer_env%astring(1:9),"(I4,A,I4)") timer_env%mepos,":",&
         timer_env%num_pe
    WRITE (timer_env%bstring(1:9),"(I4,A,I4)") timer_env%mepos,":",&
         timer_env%num_pe

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
    CHARACTER(len=max_clock_name), &
      DIMENSION(:), INTENT(in)               :: list
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

