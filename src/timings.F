!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE timings
  
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : zero
  USE stop_program, ONLY : stop_prg
  USE string_utilities, ONLY : xstring, uppercase
  USE timesl, ONLY : cputime, walltime
  USE util, ONLY : sort
  
  PRIVATE
  PUBLIC :: timeset, timestop, timeprint, trace_debug, print_stack
  
!!subroutine calling stack
  INTEGER, PARAMETER :: max_stack = 100
  INTEGER :: stack_size
  CHARACTER ( LEN = 30 ) :: routine_stack ( max_stack ), routine_name
  
!!subroutine timing and performance information
  INTEGER, PARAMETER :: max_timer = 200
  
  CHARACTER ( LEN = 30 ) :: clock_name ( max_timer )
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
 
CONTAINS

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
     " >>>>> ",TRIM(name)," started ..."

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
     CALL stop_prg ( 'timestop', 'handle > num_timer' )
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
     " <<<<< ",TRIM(routine_stack(stack_size))," ... ended"
     
!..remove name from subroutine stack
     stack_size = stack_size - 1
  END IF
  
END SUBROUTINE timestop

!******************************************************************************

SUBROUTINE timeprint ( iw, global_print_level )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER :: iw, global_print_level
  
! Locals
  INTEGER :: index ( max_timer ), i, j, ia, ie
  REAL ( dbl ) :: maxtime, mintime, perf
  CHARACTER ( LEN = 17 ) :: sname
  CHARACTER ( LEN = 6 ) :: punit
  CHARACTER ( LEN = 1 ) :: blank = ' ', ct ( 0 : 1 ) = (/ 'I', 'E' /)
  
!------------------------------------------------------------------------------
  
  IF ( num_timer > 0 .AND. global_print_level >= 0 ) THEN
     
     CALL sort ( clock_accu, num_timer, index )
     
     maxtime = clock_accu ( index ( num_timer ) )
     
     IF ( global_print_level > 10 ) THEN
        mintime = zero
     ELSE IF ( global_print_level > 5 ) THEN
        mintime = maxtime * 0.002_dbl
     ELSE
        mintime = maxtime * 0.02_dbl
     END IF
     
     WRITE ( iw, '( /, 1X, 79("-") )' )
     WRITE ( iw, '( " -", 77X, "-" )' )
     WRITE ( iw, '( " -", 32X, A, 32X, "-" )' ) ' T I M I N G '
     WRITE ( iw, '( " -", 77X, "-" )' )
     WRITE ( iw, '( 1X, 79("-"), / )' )
     WRITE ( iw, '( A, A, A )' ) ' SUBROUTINE', '        CALLS      ', &
          ' CPU TIME T        PERFORMANCE        ELAPSED TIME'
     
     DO i = num_timer, 1, -1
        
        j = index ( i )
        IF ( clock_accu ( j ) > mintime ) THEN
           
           IF ( clock_accu ( i ) > zero ) THEN
              perf = perf_count ( j ) / clock_accu ( i )
           ELSE
              perf = zero
           END IF
           
           sname = blank
           ia = 0
           ie = 0
           CALL xstring ( clock_name ( j ), ia, ie )
           ie = min ( ia+16-1, ie )
           sname ( 1 : ie-ia+1 ) = clock_name ( j ) ( ia : ie )
           punit = blank
           
           CALL xstring ( perf_unit ( j ), ia, ie )
           
           IF ( ia > ie ) THEN
              punit = '---   '
              WRITE ( iw, &
                   '( 1X, A16, I7, F15.2, 1X, A1, 12X, 1X, A6, 5X, F15.2 )' ) &
                   sname, routine_calls ( j ), clock_accu ( i ), &
                   ct ( clock_type ( j ) ), punit, wallclock_accu ( j )
           ELSE
              ie = MIN ( ia+6-1, ie )
              punit ( 6-ie+ia : 6 ) = perf_unit ( j ) ( ia : ie )
              WRITE ( iw, &
                   '( 1X, A16, I7, F15.2, 1X, A1, F12.2, 1X, A6, 5X, F15.2 )' &
                   ) &
                   sname, routine_calls ( j ), clock_accu ( i ), &
                   ct ( clock_type ( j ) ), perf, punit, wallclock_accu ( j )
           END IF
        END IF
     END DO
     
     WRITE ( iw, '( 1X, 79("-"), / )' )
  END IF
  
END SUBROUTINE timeprint

!******************************************************************************

SUBROUTINE trace_debug ( set, iw )

  IMPLICIT NONE
  
! Arguments
  CHARACTER (len=*) :: set
  INTEGER, OPTIONAL :: iw

! Locals
  INTEGER :: ia, ie
  
  CALL uppercase(set)
  CALL xstring ( set, ia, ie )
  IF ( set(ia:ie) == "START") THEN
    trace = .true.
  ELSEIF ( set(ia:ie) == "STOP" ) THEN
    trace = .false.
  ELSE
    call stop_prg ("trace_debug","unknown action")
  ENDIF
  IF ( PRESENT ( iw ) ) iod = iw

END SUBROUTINE trace_debug

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
!******************************************************************************

END MODULE timings
