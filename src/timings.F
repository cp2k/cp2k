!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002, 2004  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/timings [1.0] *
!!
!!   NAME
!!     timings
!!
!!   FUNCTION
!!     Timing routines for accounting
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     02.2004 made a stacked version (of stacks...) [Joost VandeVondele]
!!     11.2004 storable timer_envs (for f77 interface) [fawzi]
!!
!!   SOURCE
!******************************************************************************

MODULE timings

  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_cputime,&
                                             m_walltime,&
                                             print_memory
  USE mathconstants,                   ONLY: zero
  USE message_passing,                 ONLY: mp_abort,&
                                             mp_environ
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

  TYPE timer_env_type
      !PRIVATE
      INTEGER :: ref_count, id_nr
    ! subroutine calling stack
      INTEGER :: stack_size
      CHARACTER ( LEN = 80 ) :: routine_stack ( max_stack ), routine_name
      INTEGER :: handle_stack (max_stack)
    ! subroutine timing and performance information
      CHARACTER ( LEN = 80 ) :: clock_name ( max_timer )
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
      LOGICAL :: trace
   END TYPE timer_env_type

   TYPE timer_env_p_type
      TYPE(timer_env_type), POINTER :: timer_env
   END TYPE timer_env_p_type

   INTEGER, SAVE, PRIVATE :: stack_pointer = 0       ! position in timer_env_stack
   INTEGER, SAVE, PRIVATE :: last_timer_env_id = 0
   TYPE(timer_env_p_type), DIMENSION(max_timer_env_stack) :: timer_env_stack 
   CHARACTER(len=*), PARAMETER :: moduleN="timings"

!!*****
!******************************************************************************
CONTAINS

!!****f* timings/add_timer_env *
!!
!!   NAME
!!     add_timer_env
!!   FUNCTION
!!     adds the given timer_env to the top of the stack
!!   NOTES
!!     for each init_timer_env there should be the symmetric call to 
!!     rm_timer_env
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     02.2004 created [Joost VandeVondele]
!!
!!*** **********************************************************************
  SUBROUTINE add_timer_env(timer_env)
    TYPE(timer_env_type), OPTIONAL, POINTER  :: timer_env

    stack_pointer=stack_pointer+1
    IF (stack_pointer>max_timer_env_stack) THEN
       PRINT*,'stack_pointer too large is add_timer_env'
       CALL mp_abort()
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

!!****f* timings/timer_env_create *
!!
!!   NAME
!!     timer_env_create
!!   FUNCTION
!!     creates a new timer env
!!   NOTES
!!     -
!!   INPUTS
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE timer_env_create(timer_env)
    TYPE(timer_env_type), POINTER            :: timer_env

    INTEGER                                  :: stat

    ALLOCATE(timer_env,stat=stat)
    IF (stat/=0) THEN
      PRINT*,moduleN//':timer_env_create error allocationg timer_env'
      CALL mp_abort()
    END IF
    last_timer_env_id=last_timer_env_id+1
    timer_env%id_nr=last_timer_env_id
    timer_env%ref_count=1
    timer_env%stack_size=0
    timer_env%handle_stack=0
    timer_env%num_timer=0
    timer_env%trace=.FALSE.
END SUBROUTINE

!!****f* timings/rm_timer_env *
!!
!!   NAME
!!     rm_timer_env
!!   FUNCTION
!!     removes the current timer env from the stack
!!   NOTES
!!     for each rm_timer_env there should have been the symmetric call to 
!!     add_timer_env
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     02.2004 created [Joost VandeVondele]
!!
!!*** **********************************************************************
SUBROUTINE rm_timer_env()
  IF (stack_pointer<1) THEN
     PRINT*,'no timer env available in destroy_timer_env'
     CALL mp_abort()
  ENDIF
  CALL timer_env_release(timer_env_stack(stack_pointer)%timer_env)
  stack_pointer=stack_pointer-1
END SUBROUTINE

!!****f* timings/get_timer_env *
!!
!!   NAME
!!     get_timer_env
!!   FUNCTION
!!     returns the current timer env from the stack
!!   NOTES
!!     -
!!   INPUTS
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION get_timer_env() RESULT(res)
    TYPE(timer_env_type), POINTER            :: res

  IF (stack_pointer<1) THEN
     PRINT*,'no timer env available in get_timer_env'
     CALL mp_abort()
  ENDIF
  res => timer_env_stack(stack_pointer)%timer_env
END FUNCTION get_timer_env

!!****f* timings/timer_env_retain *
!!
!!   NAME
!!     timer_env_retain
!!   FUNCTION
!!     retains the given timer env
!!   NOTES
!!     -
!!   INPUTS
!!     - timer_env: the timer env to retain
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE timer_env_retain(timer_env)
    TYPE(timer_env_type), POINTER            :: timer_env

  IF (.NOT.ASSOCIATED(timer_env)) THEN
     PRINT*,moduleN,':timer_env_retain ERROR, unassociated timer_env'
     CALL mp_abort()
  ENDIF
  IF (timer_env%ref_count<1) THEN
     PRINT*,moduleN,':timer_env_retain ERROR, timer_env%ref_count=',&
          timer_env%ref_count
     CALL mp_abort()
  ENDIF
  timer_env%ref_count=timer_env%ref_count+1
END SUBROUTINE timer_env_retain

!!****f* timings/timer_env_release *
!!
!!   NAME
!!     timer_env_release
!!   FUNCTION
!!     releases the given timer env
!!   NOTES
!!     -
!!   INPUTS
!!     - timer_env: the timer env to release
!!   AUTHOR
!!     fawzi
!!
!!   SOURCE
!!*** **********************************************************************
SUBROUTINE timer_env_release(timer_env)
    TYPE(timer_env_type), POINTER            :: timer_env

    INTEGER                                  :: stat

  IF (ASSOCIATED(timer_env)) THEN
     IF (timer_env%ref_count<1) THEN
        PRINT*,moduleN,':timer_env_release ERROR, timer_env%ref_count=',&
             timer_env%ref_count,"id_nr=",timer_env%id_nr
        CALL mp_abort()
     ENDIF
     timer_env%ref_count=timer_env%ref_count-1
     IF (timer_env%ref_count==0) THEN
        timer_env%stack_size=0
        timer_env%handle_stack=0
        timer_env%num_timer=0
        timer_env%trace=.FALSE.
        DEALLOCATE(timer_env,stat=stat)
        IF (stat/=0) THEN
           PRINT*,moduleN,':timer_env_release ERROR deallocating, stat=',stat
           CALL mp_abort()
        ENDIF
     END IF
  END IF
  NULLIFY(timer_env)
END SUBROUTINE timer_env_release

!******************************************************************************
!!****** timings/timeset [1.0] *
!!
!!   NAME
!!     timeset
!!
!!   FUNCTION
!!     Start timer
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE timeset_short(name, handle ) ! provide the usual shortcut
    CHARACTER(LEN=*), INTENT(IN)             :: name
    INTEGER                                  :: handle

   CALL timeset_full(name, "I", "", handle)
END SUBROUTINE timeset_short

SUBROUTINE timeset_full ( name, ct, pu, handle )

  IMPLICIT NONE

! Arguments
  CHARACTER ( LEN = * ), INTENT ( IN ) :: name
  CHARACTER ( LEN = 1 ), INTENT ( IN ) :: ct    ! Clock type [0|1]
  CHARACTER ( LEN = * ), INTENT ( IN ) :: pu    ! Performance unit
  INTEGER, INTENT ( OUT ) :: handle

! Locals
  INTEGER :: i,stack_size
  CHARACTER ( LEN = 40 )  :: sformat
  CHARACTER ( LEN = 20 )  :: mytag
  CHARACTER ( LEN = 400 ) :: mystring
  TYPE(timer_env_type), POINTER :: timer_env
#if defined(__HPM)
#include "f_hpm.h"
#endif

!------------------------------------------------------------------------------

   IF (stack_pointer<1) THEN
      PRINT*,'no timer env available in timeset'
      CALL mp_abort()
   ENDIF

   timer_env => timer_env_stack(stack_pointer)%timer_env
!..add name to subroutine stack
  timer_env%stack_size = timer_env%stack_size + 1
  stack_size = timer_env%stack_size
  timer_env%routine_stack ( stack_size ) = name


!..is this a new timer or has this routine been called before
  handle = 0
  timer_env%routine_name = name
  DO i = 1,  timer_env%num_timer
     IF ( timer_env%routine_name == &
          timer_env%clock_name ( i ) ) handle = i
  END DO

  IF ( handle == 0 ) THEN
!..new timer started
     timer_env%num_timer = timer_env%num_timer + 1
     handle = timer_env%num_timer
     timer_env%clock_name ( handle ) = &
        timer_env%routine_name
     timer_env%perf_unit ( handle ) = pu
     IF ( ct == 'I' ) timer_env%clock_type ( handle ) = 0
     IF ( ct == 'E' ) timer_env%clock_type ( handle ) = 1
     timer_env%clock_status ( handle ) = 0
     timer_env%clock_accu ( handle ) = zero
     timer_env%wallclock_accu ( handle ) = zero
     timer_env%perf_count ( handle ) = zero
     timer_env%routine_calls ( handle ) = 0
     timer_env%routine_stack_depth ( handle ) = 0
  END IF

#if defined(__HPM)
  IF (use_HPM) THEN
      mytag=TRIM(timer_env%routine_name)
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
  WRITE(sformat,*) "(A,",MAX(1,3*stack_size-9),"X,I4,1X,I6,1X,A,A)"
  IF ( timer_env%trace ) WRITE ( mystring, sformat) &
     timer_env%astring,stack_size, &
     timer_env%routine_calls ( handle ), &
     TRIM( timer_env%routine_stack ( stack_size )),"       start"
  IF ( timer_env%trace ) CALL print_memory(TRIM(mystring))

END SUBROUTINE timeset_full

!!*****
!******************************************************************************
!!****** timings/timestop [1.0] *
!!
!!   NAME
!!     timeset
!!
!!   FUNCTION
!!     End timer
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE timestop_short( handle ) ! provide the usual shortcut
    INTEGER                                  :: handle

   CALL timestop_full(0.0_dp,handle)
END SUBROUTINE timestop_short

SUBROUTINE timestop_full ( perf, handle )

  IMPLICIT NONE

! Arguments
  REAL (KIND=dp) :: perf
  INTEGER :: handle

! Locals
  REAL (KIND=dp) :: cpt, wt
  INTEGER :: i
  CHARACTER ( LEN = 40 )  :: sformat
  CHARACTER ( LEN = 400 )  :: mystring
  TYPE(timer_env_type), POINTER :: timer_env
#if defined(__HPM)
#include "f_hpm.h"
#endif

!------------------------------------------------------------------------------

   IF (stack_pointer<1) THEN
      PRINT*,'no timer env available in timestop'
      CALL mp_abort()
   ENDIF

   timer_env => timer_env_stack(stack_pointer)%timer_env

  IF ( handle /= timer_env%handle_stack (timer_env%stack_size) ) THEN
      PRINT*,'mismatched timestop '&
        //TRIM(timer_env%routine_stack(timer_env%stack_size))//&
             ' in routine timestop'
      CALL mp_abort()
  END IF
  timer_env%handle_stack (timer_env%stack_size)=0

  IF ( handle > timer_env%num_timer ) THEN
!    CALL stop_program ( 'timestop', 'handle > num_timer' )
! intermediate work-around
     PRINT*,'handle > num_timer in routine timestop'
     CALL mp_abort()
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
  WRITE(sformat,*) "(A,",MAX(1,3*timer_env%stack_size-9),"X,I4,1X,I6,1X,A,F12.3)"
  IF ( timer_env%trace ) WRITE ( mystring, sformat) &
     timer_env%bstring, & 
     timer_env%stack_size, &
     timer_env%routine_calls ( handle ), &
     TRIM(timer_env%routine_stack(timer_env%stack_size)),cpt
  IF ( timer_env%trace ) CALL print_memory(TRIM(mystring ))

!..remove name from subroutine stack
     timer_env%stack_size = timer_env%stack_size - 1
  END IF

END SUBROUTINE timestop_full

!!*****
!******************************************************************************
!!****** timings/timeprint [1.0] *
!!
!!   NAME
!!     timeprint
!!
!!   FUNCTION
!!     Print accumulated information on timers
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE timeprint(iw,global_print_level)


    INTEGER, INTENT(IN)                      :: iw, global_print_level

    INTEGER                                  :: i, j
    CHARACTER(LEN=1), DIMENSION(0:1)         :: ct = (/"I","E"/)
    INTEGER, DIMENSION(max_timer)            :: index
    REAL(KIND=dp)                            :: asd, maxtime, mintime, perf
    TYPE(timer_env_type), POINTER            :: timer_env

!------------------------------------------------------------------------------

  IF (stack_pointer<1) THEN
      PRINT*,'no timer env available in time_set'
      CALL mp_abort()
  ENDIF

   timer_env => timer_env_stack(stack_pointer)%timer_env

  IF ((timer_env%num_timer > 0).AND.(global_print_level >= 0)) THEN

     CALL sort(timer_env%wallclock_accu,timer_env%num_timer,index)

     maxtime = timer_env%wallclock_accu(timer_env%num_timer)

     IF (global_print_level > 10) THEN
        mintime = zero
     ELSE IF (global_print_level > 5) THEN
        mintime = maxtime*0.002_dp
     ELSE
        mintime = maxtime*0.02_dp
     END IF

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

!!*****
!******************************************************************************
!!****** timings/trace_debug [1.0] *
!!
!!   NAME
!!     trace_debug
!!
!!   FUNCTION
!!     Set routine tracer
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE trace_debug ( set, groupid )


    CHARACTER(len=*)                         :: set
    INTEGER                                  :: groupid

    CHARACTER(len=5)                         :: set2
    INTEGER                                  :: numtask, taskid
    TYPE(timer_env_type), POINTER            :: timer_env

! this is needed on SUN if the trace_debug("start") is used in the call (i.e. explicit string)

   IF (stack_pointer<1) THEN
      PRINT*,'no timer env available in trace_debug'
      CALL mp_abort()
   ENDIF
 
   timer_env => timer_env_stack(stack_pointer)%timer_env

  set2=TRIM(set)
  CALL uppercase(set2)

  SELECT CASE (set2)
  CASE ("START")
    timer_env%trace = .TRUE.
  CASE ("STOP")
    timer_env%trace = .FALSE.
  CASE DEFAULT
!   CALL stop_program("trace_debug","unknown action")
!   intermediate work-around
    PRINT*,'unknown action requested in routine trace_debug'
    CALL mp_abort()
  END SELECT

  timer_env%astring = " >>>>>>>>>>>>> "
  timer_env%bstring = " >>>>>>>>>>>>> "
#if defined (__parallel)
  CALL mp_environ ( numtask, taskid, groupid )
  WRITE (timer_env%astring(1:9),"(I4,A,I4)") taskid,":",&
                  numtask
  WRITE (timer_env%bstring(1:9),"(I4,A,I4)") taskid,":",&
                  numtask
#endif

END SUBROUTINE trace_debug

!!*****
!******************************************************************************
!!****** timings/print_stack [1.0] *
!!
!!   NAME
!!     print_stack
!!
!!   FUNCTION
!!     Print current routine stack
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE print_stack ( iw )

    INTEGER                                  :: iw

    INTEGER                                  :: i, io
    TYPE(timer_env_type), POINTER            :: timer_env

  io = iw

  IF (stack_pointer<1) THEN
      PRINT*,'no timer env available in print_stack'
      CALL mp_abort()
  ENDIF

  timer_env => timer_env_stack(stack_pointer)%timer_env

  WRITE ( io , '(/,A,/)') " ===== Routine Calling Stack ===== "
  DO i = timer_env%stack_size, 1, -1
     WRITE ( io, '(T10,I4,1X,A)') i,TRIM(timer_env%routine_stack(i))
  END DO

END SUBROUTINE print_stack

! *****************************************************************************

END MODULE timings

!******************************************************************************
