!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
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
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE timings

  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE machine,                         ONLY: m_cputime,&
                                             m_flush,&
                                             m_getpid,&
                                             m_walltime,&
                                             print_memory
  USE mathconstants,                   ONLY: zero
  USE message_passing,                 ONLY: mp_abort,&
                                             mp_sync,&
                                             mp_world_setup
  USE string_utilities,                ONLY: compress,&
                                             uppercase,&
                                             xstring
  USE util,                            ONLY: sort

  PRIVATE

  PUBLIC :: timeset, timestop, timeprint, trace_debug, print_stack,&
            write_cputime, use_HPM

! want to use HPM if available ? needs care switching true/false (i.e. use_HPM=.true. timeset ... timestop use_HPM=false)
  LOGICAL :: use_HPM = .FALSE.
! subroutine calling stack
  INTEGER, PARAMETER :: max_stack = 5000
  INTEGER :: stack_size = 0
  CHARACTER ( LEN = 80 ) :: routine_stack ( max_stack ), routine_name
  INTEGER :: handle_stack (max_stack)=0

! subroutine timing and performance information
  INTEGER, PARAMETER :: max_timer = 5000

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
  INTEGER :: num_timer = 0

  INTEGER :: iod = 6
  LOGICAL :: trace = .FALSE.
  CHARACTER ( LEN = 15 ) :: astring, bstring

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** timings/timeset [1.0] *
!!
!!   NAME
!!     timeset
!!
!!   SYNOPSIS
!!     Subroutine timeset(name, ct, pu, handle)
!!       Implicit None
!!       Character(Len=*), Intent (IN):: name
!!       Character(Len=1), Intent (IN):: ct
!!       Character(Len=*), Intent (IN):: pu
!!       Integer, Intent (OUT):: handle
!!     End Subroutine timeset
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
  INTEGER :: i
  CHARACTER ( LEN = 40 )  :: sformat
  CHARACTER ( LEN = 20 )  :: mytag
  CHARACTER ( LEN = 400 ) :: mystring
#if defined(__HPM)
#include "f_hpm.h"
#endif

!------------------------------------------------------------------------------

!..add name to subroutine stack
  stack_size = stack_size + 1
  routine_stack ( stack_size ) = name


!..is this a new timer or has this routine been called before
  handle = 0
  routine_name = name
  DO i = 1, num_timer
     IF ( routine_name == clock_name ( i ) ) handle = i
  END DO

  IF ( handle == 0 ) THEN
!..new timer started
     num_timer = num_timer + 1
     handle = num_timer
     clock_name ( handle ) = routine_name
     perf_unit ( handle ) = pu
     IF ( ct == 'I' ) clock_type ( handle ) = 0
     IF ( ct == 'E' ) clock_type ( handle ) = 1
     clock_status ( handle ) = 0
     clock_accu ( handle ) = zero
     wallclock_accu ( handle ) = zero
     perf_count ( handle ) = zero
     routine_calls ( handle ) = 0
     routine_stack_depth ( handle ) = 0
  END IF

#if defined(__HPM)
  IF (use_HPM) THEN
      mytag=TRIM(routine_name)
      mytag(20:20)=CHAR(0)
      CALL f_hpmstart(handle,mytag)
  ENDIF
#endif
 
!..start clock (if not yet started, i.e. recursive calls)
  IF ( clock_status ( handle ) ==0 ) THEN
     clock_start ( handle ) = m_cputime ( )
     wallclock_start ( handle ) = m_walltime ( )
     clock_status ( handle ) = 1
  END IF

  routine_calls ( handle ) = routine_calls ( handle ) + 1
  routine_stack_depth ( handle ) = routine_stack_depth ( handle ) + stack_size

  handle_stack (stack_size) = handle

!..if debug mode echo the subroutine name
  WRITE(sformat,*) "(A,",MAX(1,3*stack_size-9),"X,I4,1X,I6,1X,A,A)"
  IF ( trace ) WRITE ( mystring, sformat) &
     astring,stack_size,routine_calls ( handle ),TRIM( routine_stack ( stack_size )),"       start"
  IF ( trace ) CALL print_memory(TRIM(mystring))

END SUBROUTINE timeset

!!*****
!******************************************************************************
!!****** timings/timeset [1.0] *
!!
!!   NAME
!!     timeset
!!
!!   SYNOPSIS
!!     Subroutine timeset(name, ct, pu, handle)
!!       Implicit None
!!       Character(Len=*), Intent (IN):: name
!!       Character(Len=1), Intent (IN):: ct
!!       Character(Len=*), Intent (IN):: pu
!!       Integer, Intent (OUT):: handle
!!     End Subroutine timeset
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

  IF ( handle /= handle_stack (stack_size) ) THEN
!     CALL stop_program ( 'timestop', 'mismatched timestop '//&
!          TRIM(routine_stack(stack_size)) )
! intermediate work-around
      print*,'mismatched timestop '//TRIM(routine_stack(stack_size))//&
             ' in routine timestop'
      CALL mp_abort()
  END IF
  handle_stack (stack_size)=0

  IF ( handle > num_timer ) THEN
!    CALL stop_program ( 'timestop', 'handle > num_timer' )
! intermediate work-around
     print*,'handle > num_timer in routine timestop'
     CALL mp_abort()
  ELSE
#if defined(__HPM)
     IF (use_HPM) THEN
         CALL f_hpmstop(handle)
     ENDIF
#endif
     cpt = m_cputime ( ) - clock_start ( handle )
     wt = m_walltime ( ) - wallclock_start ( handle )
     clock_accu ( handle ) = clock_accu ( handle ) + cpt
     wallclock_accu ( handle ) = wallclock_accu ( handle ) + wt
     clock_status ( handle ) = 0
     perf_count ( handle ) = perf_count ( handle ) + perf

     DO i = 1, num_timer
        IF ( clock_status ( i ) == 1 .AND. clock_type ( i ) == 1 ) THEN
           clock_start ( i ) = clock_start ( i ) + cpt
           wallclock_start ( i ) = wallclock_start ( i ) + wt
        END IF
     END DO

!..if debug mode echo the subroutine name
  WRITE(sformat,*) "(A,",MAX(1,3*stack_size-9),"X,I4,1X,I6,1X,A,F12.3)"
  IF ( trace ) WRITE ( mystring, sformat) &
     bstring,stack_size,routine_calls ( handle ),TRIM(routine_stack(stack_size)),cpt
  IF ( trace ) CALL print_memory(TRIM(mystring ))

!..remove name from subroutine stack
     stack_size = stack_size - 1
  END IF

END SUBROUTINE timestop

!!*****
!******************************************************************************
!!****** timings/timeprint [1.0] *
!!
!!   NAME
!!     timeprint
!!
!!   SYNOPSIS
!!     Subroutine timeprint(iw, global_print_level)
!!       Implicit None
!!       Integer, Intent (IN):: global_print_level, iw
!!     End Subroutine timeprint
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

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT(IN) :: global_print_level,iw

! Locals
  REAL(dbl) :: maxtime,mintime,perf,asd
  INTEGER   :: i,j

  CHARACTER(LEN=1), DIMENSION(0:1) :: ct = (/"I","E"/)

  INTEGER, DIMENSION(max_timer) :: index

!------------------------------------------------------------------------------

  IF ((num_timer > 0).AND.(global_print_level >= 0)) THEN

     CALL sort(wallclock_accu,num_timer,index)

     maxtime = wallclock_accu(num_timer)

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

     DO i=num_timer,1,-1

       j = INDEX(i)

       IF (wallclock_accu(i) > mintime) THEN

         IF (wallclock_accu(i) > zero) THEN
           perf = perf_count(j)/wallclock_accu(i)
         ELSE
           perf = zero
         END IF
         ! average stack depth
         asd=DBLE(routine_stack_depth(j))/DBLE(MAX(1,routine_calls(j)))

         IF (LEN_TRIM(perf_unit(j)) == 0) THEN
           WRITE (UNIT=iw,&
                  FMT="(T2,A,T30,I7,T38,F4.1,T42,F9.2,T52,A,T60,A,T68,F13.2)")&
             TRIM(clock_name(j)),routine_calls(j),ASD,clock_accu(j),&
             ct(clock_type(j)),"---",wallclock_accu(i)
         ELSE
           WRITE (UNIT=iw,&
                  FMT="(T2,A,T30,I7,T38,F4.1,T42,F9.2,T52,A,T54,F5.1,T60,A,T68,F13.2)")&
             TRIM(clock_name(j)),routine_calls(j),ASD,clock_accu(j),&
             ct(clock_type(j)),perf,TRIM(perf_unit(j)),wallclock_accu(i)
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
!!   SYNOPSIS
!!     Subroutine trace_debug(set, iw)
!!       Implicit None
!!       Character(Len=*):: set
!!       Integer, Optional:: iw
!!     End Subroutine trace_debug
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

SUBROUTINE trace_debug ( set, iw )

  IMPLICIT NONE

! Arguments 
  CHARACTER (len=*) :: set
  CHARACTER (len=5) :: set2
  INTEGER, OPTIONAL :: iw

! Locals
  INTEGER :: numtask, taskid, groupid
  ! this is needed on SUN if the trace_debug("start") is used in the call (i.e. explicit string)
  set2=TRIM(set)
  CALL uppercase(set2)

  SELECT CASE (set2)
  CASE ("START")
    trace = .TRUE.
  CASE ("STOP")
    trace = .FALSE.
  CASE DEFAULT
!   CALL stop_program("trace_debug","unknown action")
!   intermediate work-around
    print*,'unknown action requested in routine trace_debug'
    CALL mp_abort()
  END SELECT

  IF ( PRESENT ( iw ) ) iod = iw

  astring = " >>>>>>>>>>>>> "
  bstring = " >>>>>>>>>>>>> "
#if defined (__parallel)
  CALL mp_world_setup ( numtask, taskid, groupid )
  WRITE (astring(1:9),"(I4,A,I4)") taskid,":",numtask
  WRITE (bstring(1:9),"(I4,A,I4)") taskid,":",numtask
#endif

END SUBROUTINE trace_debug

!!*****
!******************************************************************************
!!****** timings/print_stack [1.0] *
!!
!!   NAME
!!     print_stack
!!
!!   SYNOPSIS
!!     Subroutine print_stack(iw)
!!       Implicit None
!!       Integer, Optional:: iw
!!     End Subroutine print_stack
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

  IMPLICIT NONE

! Arguments
  INTEGER, OPTIONAL :: iw

! Locals
  INTEGER :: io, i

  io = 0
  IF ( PRESENT ( iw ) ) io = iw

  WRITE ( io , '(/,A,/)') " ===== Routine Calling Stack ===== "
  DO i = stack_size, 1, -1
     WRITE ( io, '(T10,I4,1X,A)') i,TRIM(routine_stack(i))
  END DO

END SUBROUTINE print_stack

! *****************************************************************************

  SUBROUTINE write_cputime(checkpoint,globenv)

!   Purpose: Write the current CPU time of all processes in seconds.

!   History: - Creation (25.09.2002,MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    CHARACTER(LEN=*), INTENT(IN)              :: checkpoint

!   *** Local variables ***

    CHARACTER(LEN=LEN_TRIM(checkpoint)+7)  :: string1
    CHARACTER(LEN=LEN_TRIM(checkpoint)+55) :: string2
    INTEGER                                :: group,ipe,mype,npe,output_unit,&
                                              pid

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%print%cputime) RETURN

    group = globenv%group
    npe = globenv%num_pe
    mype = globenv%mepos
    output_unit = globenv%scr

    CALL m_getpid(pid)

    DO ipe=0,npe-1
      IF (ipe == mype) THEN
        WRITE (UNIT=string1,FMT="(I6,A)") ipe,":"//TRIM(checkpoint)
        string1 = ADJUSTL(string1)
        WRITE (UNIT=string2,FMT="(A,I8,A,F12.3,A)")&
          "CPU time of process",pid," <"//TRIM(string1)//">: ",m_cputime(),&
          " sec"
        CALL compress(string2)
        WRITE (UNIT=output_unit,FMT="(/,T2,A)",ADVANCE="NO") TRIM(string2)
      END IF
      CALL m_flush(output_unit)
      CALL mp_sync(group)
    END DO

  END SUBROUTINE write_cputime

!!*****
!******************************************************************************

END MODULE timings

!******************************************************************************
