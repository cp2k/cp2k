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
!!
!!   SOURCE
!******************************************************************************

MODULE timings

  USE kinds,                           ONLY: dbl
  USE machine,                         ONLY: m_cputime,&
                                             m_walltime,&
                                             print_memory
  USE mathconstants,                   ONLY: zero
  USE message_passing,                 ONLY: mp_abort,&
                                             mp_environ
  USE string_utilities,                ONLY: uppercase
  USE util,                            ONLY: sort

  PRIVATE

  PUBLIC :: timeset, timestop, timeprint, trace_debug, print_stack,&
            use_HPM

  ! used to create a fresh timer_env
  PUBLIC :: init_timer_env, destroy_timer_env

! want to use HPM if available ? needs care switching true/false (i.e. use_HPM=.true. timeset ... timestop use_HPM=false)
  LOGICAL :: use_HPM = .FALSE.

! hard coded stack sizes
  INTEGER, PARAMETER :: max_timer_env_stack=10  ! number of stacked timers
  INTEGER, PARAMETER :: max_stack = 50          ! the calling stack max depth for a given timer env
  INTEGER, PARAMETER :: max_timer = 2000        ! total number of different callers

  TYPE timer_env_type
    ! subroutine calling stack
      INTEGER :: stack_size
      CHARACTER ( LEN = 80 ) :: routine_stack ( max_stack ), routine_name
      INTEGER :: handle_stack (max_stack)
    ! subroutine timing and performance information
      CHARACTER ( LEN = 80 ) :: clock_name ( max_timer )
      INTEGER :: routine_calls ( max_timer )
      INTEGER :: routine_stack_depth ( max_timer )
      INTEGER :: clock_status ( max_timer )
      REAL ( dbl ) :: clock_start ( max_timer )
      REAL ( dbl ) :: clock_accu ( max_timer )
      REAL ( dbl ) :: wallclock_start ( max_timer )
      REAL ( dbl ) :: wallclock_accu ( max_timer )
      INTEGER :: clock_type ( max_timer )
      REAL ( dbl ) :: perf_count ( max_timer )
      CHARACTER ( LEN = 10 ) :: perf_unit ( max_timer )
      INTEGER :: num_timer 
      CHARACTER ( LEN = 15 ) :: astring, bstring
      LOGICAL :: trace
   END TYPE timer_env_type

   INTEGER              :: stack_pointer = 0       ! position in timer_env_stack
   TYPE(timer_env_type), DIMENSION(max_timer_env_stack) :: timer_env_stack 

!!*****
!******************************************************************************
CONTAINS

!!****f* timings/init_timer_env *
!!
!!   NAME
!!     init_timer_env
!!   FUNCTION
!!     creates a fresh timer env (on the stack)
!!   NOTES
!!     for each init_timer_env there should be the symmetric call to destroy timer env
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     02.2004 created [Joost VandeVondele]
!!
!!   SOURCE
!!*** **********************************************************************
SUBROUTINE init_timer_env()
    stack_pointer=stack_pointer+1
    IF (stack_pointer>max_timer_env_stack) THEN
      PRINT*,'stack_pointer too large is init_timer_env'
      CALL mp_abort()
    ENDIF
    timer_env_stack(stack_pointer)%stack_size=0
    timer_env_stack(stack_pointer)%handle_stack=0
    timer_env_stack(stack_pointer)%num_timer=0
    timer_env_stack(stack_pointer)%trace=.FALSE.
END SUBROUTINE
!!****f* timings/destroy_timer_env *
!!
!!   NAME
!!     destroy_timer_env
!!   FUNCTION
!!     invalidates the current timer env
!!   NOTES
!!     for each init_timer_env there should be the symmetric call to destroy timer env
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     02.2004 created [Joost VandeVondele]
!!
!!   SOURCE
!!*** **********************************************************************
SUBROUTINE destroy_timer_env()
    timer_env_stack(stack_pointer)%stack_size=0
    timer_env_stack(stack_pointer)%handle_stack=0
    timer_env_stack(stack_pointer)%num_timer=0
    timer_env_stack(stack_pointer)%trace=.FALSE.
    stack_pointer=stack_pointer-1
    IF (stack_pointer<0) THEN
      PRINT*,'no timer env available in destroy_timer_env'
      CALL mp_abort()
    ENDIF
END SUBROUTINE


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

SUBROUTINE timeset ( name, ct, pu, handle )

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
#if defined(__HPM)
#include "f_hpm.h"
#endif

!------------------------------------------------------------------------------

   IF (stack_pointer<1) THEN
      PRINT*,'no timer env available in timeset'
      CALL mp_abort()
   ENDIF


!..add name to subroutine stack
  timer_env_stack(stack_pointer)%stack_size = timer_env_stack(stack_pointer)%stack_size + 1
  stack_size = timer_env_stack(stack_pointer)%stack_size
  timer_env_stack(stack_pointer)%routine_stack ( stack_size ) = name


!..is this a new timer or has this routine been called before
  handle = 0
  timer_env_stack(stack_pointer)%routine_name = name
  DO i = 1,  timer_env_stack(stack_pointer)%num_timer
     IF ( timer_env_stack(stack_pointer)%routine_name == &
          timer_env_stack(stack_pointer)%clock_name ( i ) ) handle = i
  END DO

  IF ( handle == 0 ) THEN
!..new timer started
     timer_env_stack(stack_pointer)%num_timer = timer_env_stack(stack_pointer)%num_timer + 1
     handle = timer_env_stack(stack_pointer)%num_timer
     timer_env_stack(stack_pointer)%clock_name ( handle ) = &
        timer_env_stack(stack_pointer)%routine_name
     timer_env_stack(stack_pointer)%perf_unit ( handle ) = pu
     IF ( ct == 'I' ) timer_env_stack(stack_pointer)%clock_type ( handle ) = 0
     IF ( ct == 'E' ) timer_env_stack(stack_pointer)%clock_type ( handle ) = 1
     timer_env_stack(stack_pointer)%clock_status ( handle ) = 0
     timer_env_stack(stack_pointer)%clock_accu ( handle ) = zero
     timer_env_stack(stack_pointer)%wallclock_accu ( handle ) = zero
     timer_env_stack(stack_pointer)%perf_count ( handle ) = zero
     timer_env_stack(stack_pointer)%routine_calls ( handle ) = 0
     timer_env_stack(stack_pointer)%routine_stack_depth ( handle ) = 0
  END IF

#if defined(__HPM)
  IF (use_HPM) THEN
      mytag=TRIM(timer_env_stack(stack_pointer)%routine_name)
      mytag(20:20)=CHAR(0)
      CALL f_hpmstart(handle,mytag)
  ENDIF
#endif
 
!..start clock (if not yet started, i.e. recursive calls)
  IF ( timer_env_stack(stack_pointer)%clock_status ( handle ) ==0 ) THEN
     timer_env_stack(stack_pointer)%clock_start ( handle ) = m_cputime ( )
     timer_env_stack(stack_pointer)%wallclock_start ( handle ) = m_walltime ( )
     timer_env_stack(stack_pointer)%clock_status ( handle ) = 1
  END IF

  timer_env_stack(stack_pointer)%routine_calls ( handle ) = &
         timer_env_stack(stack_pointer)%routine_calls ( handle ) + 1
  timer_env_stack(stack_pointer)%routine_stack_depth ( handle ) =&
         timer_env_stack(stack_pointer)%routine_stack_depth ( handle ) + stack_size

  timer_env_stack(stack_pointer)%handle_stack (stack_size) = handle

!..if debug mode echo the subroutine name
  WRITE(sformat,*) "(A,",MAX(1,3*stack_size-9),"X,I4,1X,I6,1X,A,A)"
  IF ( timer_env_stack(stack_pointer)%trace ) WRITE ( mystring, sformat) &
     timer_env_stack(stack_pointer)%astring,stack_size, &
     timer_env_stack(stack_pointer)%routine_calls ( handle ), &
     TRIM( timer_env_stack(stack_pointer)%routine_stack ( stack_size )),"       start"
  IF ( timer_env_stack(stack_pointer)%trace ) CALL print_memory(TRIM(mystring))

END SUBROUTINE timeset

!!*****
!******************************************************************************
!!****** timings/timeset [1.0] *
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

SUBROUTINE timestop ( perf, handle )

  IMPLICIT NONE

! Arguments
  REAL ( dbl ) :: perf
  INTEGER :: handle

! Locals
  REAL ( dbl ) :: cpt, wt
  INTEGER :: i
  CHARACTER ( LEN = 40 )  :: sformat
  CHARACTER ( LEN = 400 )  :: mystring
#if defined(__HPM)
#include "f_hpm.h"
#endif

!------------------------------------------------------------------------------

   IF (stack_pointer<1) THEN
      PRINT*,'no timer env available in timestop'
      CALL mp_abort()
   ENDIF


  IF ( handle /= timer_env_stack(stack_pointer)%handle_stack (timer_env_stack(stack_pointer)%stack_size) ) THEN
      PRINT*,'mismatched timestop '&
        //TRIM(timer_env_stack(stack_pointer)%routine_stack(timer_env_stack(stack_pointer)%stack_size))//&
             ' in routine timestop'
      CALL mp_abort()
  END IF
  timer_env_stack(stack_pointer)%handle_stack (timer_env_stack(stack_pointer)%stack_size)=0

  IF ( handle > timer_env_stack(stack_pointer)%num_timer ) THEN
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
     cpt = m_cputime ( ) - timer_env_stack(stack_pointer)%clock_start ( handle )
     wt = m_walltime ( ) - timer_env_stack(stack_pointer)%wallclock_start ( handle )
     timer_env_stack(stack_pointer)%clock_accu ( handle ) = &
            timer_env_stack(stack_pointer)%clock_accu ( handle ) + cpt
     timer_env_stack(stack_pointer)%wallclock_accu ( handle ) = &
            timer_env_stack(stack_pointer)%wallclock_accu ( handle ) + wt
     timer_env_stack(stack_pointer)%clock_status ( handle ) = 0
     timer_env_stack(stack_pointer)%perf_count ( handle ) = &
            timer_env_stack(stack_pointer)%perf_count ( handle ) + perf

     DO i = 1, timer_env_stack(stack_pointer)% num_timer
        IF ( timer_env_stack(stack_pointer)%clock_status ( i ) == 1 .AND. &
             timer_env_stack(stack_pointer)%clock_type ( i ) == 1 ) THEN
           timer_env_stack(stack_pointer)%clock_start ( i ) = &
               timer_env_stack(stack_pointer)%clock_start ( i ) + cpt
           timer_env_stack(stack_pointer)%wallclock_start ( i ) = &
               timer_env_stack(stack_pointer)%wallclock_start ( i ) + wt
        END IF
     END DO

!..if debug mode echo the subroutine name
  WRITE(sformat,*) "(A,",MAX(1,3*timer_env_stack(stack_pointer)%stack_size-9),"X,I4,1X,I6,1X,A,F12.3)"
  IF ( timer_env_stack(stack_pointer)%trace ) WRITE ( mystring, sformat) &
     timer_env_stack(stack_pointer)%bstring, & 
     timer_env_stack(stack_pointer)%stack_size, &
     timer_env_stack(stack_pointer)%routine_calls ( handle ), &
     TRIM(timer_env_stack(stack_pointer)%routine_stack(timer_env_stack(stack_pointer)%stack_size)),cpt
  IF ( timer_env_stack(stack_pointer)%trace ) CALL print_memory(TRIM(mystring ))

!..remove name from subroutine stack
     timer_env_stack(stack_pointer)%stack_size = timer_env_stack(stack_pointer)%stack_size - 1
  END IF

END SUBROUTINE timestop

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
    REAL(dbl)                                :: asd, maxtime, mintime, perf

!------------------------------------------------------------------------------

  IF (stack_pointer<1) THEN
      PRINT*,'no timer env available in time_set'
      CALL mp_abort()
  ENDIF


  IF ((timer_env_stack(stack_pointer)%num_timer > 0).AND.(global_print_level >= 0)) THEN

     CALL sort(timer_env_stack(stack_pointer)%wallclock_accu,timer_env_stack(stack_pointer)%num_timer,index)

     maxtime = timer_env_stack(stack_pointer)%wallclock_accu(timer_env_stack(stack_pointer)%num_timer)

     IF (global_print_level > 10) THEN
        mintime = zero
     ELSE IF (global_print_level > 5) THEN
        mintime = maxtime*0.002_dbl
     ELSE
        mintime = maxtime*0.02_dbl
     END IF

     WRITE (UNIT=iw,FMT="(/,T2,A)") REPEAT("-",79)
     WRITE (UNIT=iw,FMT="(T2,A,T80,A)") "-","-"
     WRITE (UNIT=iw,FMT="(T2,A,T35,A,T80,A)") "-","T I M I N G","-"
     WRITE (UNIT=iw,FMT="(T2,A,T80,A)") "-","-"
     WRITE (UNIT=iw,FMT="(T2,A)") REPEAT("-",79)
     WRITE (UNIT=iw,FMT="(T2,A,T32,A,T38,A,T43,A,T55,A,T69,A)")&
       "SUBROUTINE","CALLS"," ASD","CPU TIME T","PERFORMANCE","ELAPSED TIME"

     DO i=timer_env_stack(stack_pointer)%num_timer,1,-1

       j = INDEX(i)

       IF (timer_env_stack(stack_pointer)%wallclock_accu(i) > mintime) THEN

         IF (timer_env_stack(stack_pointer)%wallclock_accu(i) > zero) THEN
           perf = timer_env_stack(stack_pointer)%perf_count(j)/timer_env_stack(stack_pointer)%wallclock_accu(i)
         ELSE
           perf = zero
         END IF
         ! average stack depth
         asd=DBLE(timer_env_stack(stack_pointer)%routine_stack_depth(j))/ &
             DBLE(MAX(1,timer_env_stack(stack_pointer)%routine_calls(j)))

         IF (LEN_TRIM(timer_env_stack(stack_pointer)%perf_unit(j)) == 0) THEN
           WRITE (UNIT=iw,&
                  FMT="(T2,A,T30,I7,T38,F4.1,T42,F9.2,T52,A,T60,A,T68,F13.2)")&
             TRIM(timer_env_stack(stack_pointer)%clock_name(j)),timer_env_stack(stack_pointer)%routine_calls(j), & 
             ASD,timer_env_stack(stack_pointer)%clock_accu(j),&
             ct(timer_env_stack(stack_pointer)%clock_type(j)),"---",timer_env_stack(stack_pointer)%wallclock_accu(i)
         ELSE
           WRITE (UNIT=iw,&
                  FMT="(T2,A,T30,I7,T38,F4.1,T42,F9.2,T52,A,T54,F5.1,T60,A,T68,F13.2)")&
             TRIM(timer_env_stack(stack_pointer)%clock_name(j)),timer_env_stack(stack_pointer)%routine_calls(j),&
             ASD,timer_env_stack(stack_pointer)%clock_accu(j),&
             ct(timer_env_stack(stack_pointer)%clock_type(j)),perf,TRIM(timer_env_stack(stack_pointer)%perf_unit(j)),&
             timer_env_stack(stack_pointer)%wallclock_accu(i)
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

! this is needed on SUN if the trace_debug("start") is used in the call (i.e. explicit string)

   IF (stack_pointer<1) THEN
      PRINT*,'no timer env available in trace_debug'
      CALL mp_abort()
   ENDIF
 
  set2=TRIM(set)
  CALL uppercase(set2)

  SELECT CASE (set2)
  CASE ("START")
    timer_env_stack(stack_pointer)%trace = .TRUE.
  CASE ("STOP")
    timer_env_stack(stack_pointer)%trace = .FALSE.
  CASE DEFAULT
!   CALL stop_program("trace_debug","unknown action")
!   intermediate work-around
    PRINT*,'unknown action requested in routine trace_debug'
    CALL mp_abort()
  END SELECT

  timer_env_stack(stack_pointer)%astring = " >>>>>>>>>>>>> "
  timer_env_stack(stack_pointer)%bstring = " >>>>>>>>>>>>> "
#if defined (__parallel)
  CALL mp_environ ( numtask, taskid, groupid )
  WRITE (timer_env_stack(stack_pointer)%astring(1:9),"(I4,A,I4)") taskid,":",&
                  numtask
  WRITE (timer_env_stack(stack_pointer)%bstring(1:9),"(I4,A,I4)") taskid,":",&
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

  io = iw

  IF (stack_pointer<1) THEN
      PRINT*,'no timer env available in print_stack'
      CALL mp_abort()
  ENDIF


  WRITE ( io , '(/,A,/)') " ===== Routine Calling Stack ===== "
  DO i = timer_env_stack(stack_pointer)%stack_size, 1, -1
     WRITE ( io, '(T10,I4,1X,A)') i,TRIM(timer_env_stack(stack_pointer)%routine_stack(i))
  END DO

END SUBROUTINE print_stack

! *****************************************************************************

END MODULE timings

!******************************************************************************
