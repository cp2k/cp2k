!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
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

  USE kinds, ONLY : dbl

  USE mathconstants, ONLY : zero
  USE string_utilities, ONLY : xstring, uppercase
  USE termination, ONLY : stop_program
  USE timesl, ONLY : cputime, walltime
  USE util, ONLY : sort
  USE message_passing, ONLY : mp_world_setup

  PRIVATE
  PUBLIC :: timeset, timestop, timeprint, trace_debug, print_stack

! subroutine calling stack
  INTEGER, PARAMETER :: max_stack = 10000
  INTEGER :: stack_size
  CHARACTER ( LEN = 31 ) :: routine_stack ( max_stack ), routine_name

! subroutine timing and performance information
  INTEGER, PARAMETER :: max_timer = 200

  CHARACTER ( LEN = 31 ) :: clock_name ( max_timer )
  INTEGER :: routine_calls ( max_timer )
  INTEGER :: clock_status ( max_timer )
  REAL ( dbl ) :: clock_start ( max_timer )
  REAL ( dbl ) :: clock_accu ( max_timer )
  REAL ( dbl ) :: wallclock_start ( max_timer )
  REAL ( dbl ) :: wallclock_accu ( max_timer )
  INTEGER :: clock_type ( max_timer )
  REAL ( dbl ) :: perf_count ( max_timer )
  CHARACTER ( LEN = 10 ) :: perf_unit ( max_timer )
  INTEGER :: num_timer = 0

  INTEGER :: iod = 0
  LOGICAL :: trace = .false.
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

!------------------------------------------------------------------------------

!..add name to subroutine stack
  stack_size = stack_size + 1
  routine_stack ( stack_size ) = name

!..if debug mode echo the subroutine name
  IF ( trace ) write ( iod, '(A,A,A)') &
     astring,TRIM(name)," started ..."

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
  END IF

!..start clock (if not yet started, i.e. recursive calls)
  IF ( clock_status ( handle ) ==0 ) THEN
     clock_start ( handle ) = cputime ( )
     wallclock_start ( handle ) = walltime ( )
     clock_status ( handle ) = 1
  END IF

  routine_calls ( handle ) = routine_calls ( handle ) + 1

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

!------------------------------------------------------------------------------

  IF ( handle > num_timer ) THEN
     CALL stop_program ( 'timestop', 'handle > num_timer' )
  ELSE
     cpt = cputime ( ) - clock_start ( handle )
     wt = walltime ( ) - wallclock_start ( handle )
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
  IF ( trace ) write ( iod, '(A,A,A)') &
     bstring,TRIM(routine_stack(stack_size))," ... ended"

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
  REAL(dbl) :: maxtime,mintime,perf
  INTEGER   :: i,j

  CHARACTER(LEN=1), DIMENSION(0:1) :: ct = (/"I","E"/)

  INTEGER, DIMENSION(max_timer) :: index

!------------------------------------------------------------------------------

  IF ((num_timer > 0).AND.(global_print_level >= 0)) THEN

     CALL sort(clock_accu,num_timer,index)

     maxtime = clock_accu(num_timer)

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
     WRITE (UNIT=iw,FMT="(T2,A,T36,A,T43,A,T55,A,T69,A)")&
       "SUBROUTINE","CALLS","CPU TIME T","PERFORMANCE","ELAPSED TIME"

     DO i=num_timer,1,-1

       j = index(i)

       IF (clock_accu(i) > mintime) THEN

         IF (clock_accu(i) > zero) THEN
           perf = perf_count(j)/clock_accu(i)
         ELSE
           perf = zero
         END IF

         IF (LEN_TRIM(perf_unit(j)) == 0) THEN
           WRITE (UNIT=iw,&
                  FMT="(T2,A,T34,I7,T42,F9.2,T52,A,T60,A,T68,F13.2)")&
             TRIM(clock_name(j)),routine_calls(j),clock_accu(i),&
             ct(clock_type(j)),"---",wallclock_accu(j)
         ELSE
           WRITE (UNIT=iw,&
                  FMT="(T2,A,T34,I7,T42,F9.2,T52,A,T54,F5.1,T60,A,T68,F13.2)")&
             TRIM(clock_name(j)),routine_calls(j),clock_accu(i),&
             ct(clock_type(j)),perf,TRIM(perf_unit(j)),wallclock_accu(j)
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

SUBROUTINE trace_debug ( set, iw )

  IMPLICIT NONE

! Arguments
  CHARACTER (len=*) :: set
  INTEGER, OPTIONAL :: iw

! Locals
  INTEGER :: ia, ie
  INTEGER :: numtask, taskid, groupid

  CALL uppercase(set)
  CALL xstring ( set, ia, ie )
  IF ( set(ia:ie) == "START") THEN
    trace = .true.
  ELSEIF ( set(ia:ie) == "STOP" ) THEN
    trace = .false.
  ELSE
    call stop_program ("trace_debug","unknown action")
  ENDIF
  IF ( PRESENT ( iw ) ) iod = iw
#if defined (__parallel)
  CALL mp_world_setup ( numtask, taskid, groupid )
  astring = " >>>>>>>>>>>>> "
  bstring = " >>>>>>>>>>>>> "
  WRITE(astring(1:9),'(i4,a,i4)') taskid,":",numtask
  WRITE(bstring(1:9),'(i4,a,i4)') taskid,":",numtask
#else
  astring = " >>>>>>>>>>>>> "
  bstring = " >>>>>>>>>>>>> "
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

  IMPLICIT NONE

! Arguments
  INTEGER, OPTIONAL :: iw

! Locals
  INTEGER :: io, i

  io = 0
  IF ( PRESENT ( iw ) ) io = iw

  WRITE ( io , '(/,A)') " ===== Routine Calling Stack ===== "
  DO i = stack_size, 1, -1
     WRITE ( io, '(T10,A)') ADJUSTL( routine_stack ( i ) )
  END DO

END SUBROUTINE print_stack

!!*****
!******************************************************************************

END MODULE timings

!******************************************************************************
