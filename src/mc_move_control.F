!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2003  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/mc_move_control [1.0] *
!!
!!   NAME
!!     qs_mc
!!
!!   FUNCTION
!!     the various moves in Monte Carlo (MC) simulations
!!
!!   AUTHOR
!!     Matthew J. McGrath  (10.16.2003)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE mc_move_control
  USE cp_error_handling,               ONLY: cp_error_type
  USE kinds,                           ONLY: dp,&
                                             dp
  USE mathconstants,                   ONLY: pi
  USE mc_types,                        ONLY: mc_moves_type,&
                                             mc_simulation_parameters_type,&
                                             get_mc_par,&
                                             set_mc_par
  USE physcon,                         ONLY: angstrom
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_move_control"

  PUBLIC :: init_mc_moves, &
           mc_move_update,move_q_reinit,q_move_accept,mc_moves_release,&
           write_move_stats
 
! *****************************************************************************

CONTAINS

! *****************************************************************************
!!****s* mc_move_control.F/init_mc_moves
!!
!!   NAME
!!     init_mc_moves
!!
!!   FUNCTION
!!     allocates and initializes the structure to record all move
!!     attempts/successes
!!
!!   ARGUMENTS
!!     - moves: the move structure to update
!!  
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
  SUBROUTINE init_mc_moves ( moves )

    TYPE(mc_moves_type), POINTER             :: moves

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "init_mc_moves"

    INTEGER                                  :: handle, stat

!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

    CALL timeset(routineN,'I','',handle)

! allocate all the structures
      ALLOCATE (moves,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves",0)
      ALLOCATE (moves%bond,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%bond",0)
      ALLOCATE (moves%angle,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%angle",0)
      ALLOCATE (moves%trans,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%trans",0)
      ALLOCATE (moves%rot,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%rot",0)
      ALLOCATE (moves%classical_bond,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%classical_bond",0)
      ALLOCATE (moves%classical_angle,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%classical_angle",0)
      ALLOCATE (moves%classical_trans,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%classical_trans",0)
      ALLOCATE (moves%classical_rot,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%classical_rot",0)
      ALLOCATE (moves%classical_volume,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%classical_volume",0)
      ALLOCATE (moves%volume,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%volume",0)
      ALLOCATE (moves%swap,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%swap",0)
      ALLOCATE (moves%Quickstep,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%Quickstep",0)
! set all the counters equal to zero
      moves%classical_bond%attempts=0
      moves%classical_bond%successes=0
      moves%classical_angle%attempts=0
      moves%classical_angle%successes=0
      moves%classical_trans%attempts=0
      moves%classical_trans%successes=0
      moves%classical_rot%attempts=0
      moves%classical_rot%successes=0
      moves%classical_volume%attempts=0
      moves%classical_volume%successes=0

      moves%bond%attempts=0
      moves%bond%successes=0
      moves%bond%qsuccesses=0
      moves%angle%attempts=0
      moves%angle%successes=0
      moves%angle%qsuccesses=0
      moves%trans%attempts=0
      moves%trans%successes=0
      moves%trans%qsuccesses=0
      moves%rot%attempts=0
      moves%rot%successes=0
      moves%rot%qsuccesses=0
      moves%volume%attempts=0
      moves%volume%successes=0
      moves%volume%qsuccesses=0
      moves%swap%attempts=0
      moves%swap%successes=0
      moves%swap%qsuccesses=0
      moves%Quickstep%attempts=0
      moves%Quickstep%successes=0
      moves%Quickstep%qsuccesses=0
      moves%classical_energy=0.0E0_dp
      moves%trans_dis=0.0E0_dp
      moves%qtrans_dis=0.0E0_dp
      moves%empty=0
      moves%grown=0
!      moves%force_create=1

 ! end the timing
      CALL timestop(0.0E0_dp,handle)

 END SUBROUTINE init_mc_moves

! *****************************************************************************
!!****s* mc_move_control.F/mc_moves_release
!!
!!   NAME
!!     mc_moves_release
!!
!!   FUNCTION
!!     deallocates all the structures and nullifies the pointer
!!
!!   ARGUMENTS
!!     - moves: the move structure to release
!!  
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
  SUBROUTINE mc_moves_release ( moves )

    TYPE(mc_moves_type), POINTER             :: moves

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "mc_moves_release"

    INTEGER                                  :: handle, stat

!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

    CALL timeset(routineN,'I','',handle)

! allocate all the structures
      DEALLOCATE (moves%bond,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%bond",0)
      DEALLOCATE (moves%angle,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%angle",0)
      DEALLOCATE (moves%trans,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%trans",0)
      DEALLOCATE (moves%rot,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%rot",0)
      DEALLOCATE (moves%classical_bond,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%classical_bond",0)
      DEALLOCATE (moves%classical_angle,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%classical_angle",0)
      DEALLOCATE (moves%classical_trans,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%classical_trans",0)
      DEALLOCATE (moves%classical_rot,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%classical_rot",0)
      DEALLOCATE (moves%classical_volume,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%classical_volume",0)
      DEALLOCATE (moves%volume,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%volume",0)
      DEALLOCATE (moves%swap,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%swap",0)
      DEALLOCATE (moves%Quickstep,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%Quickstep",0)

      DEALLOCATE (moves,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves",0)

! now nullify the moves
      NULLIFY(moves)
      
 ! end the timing
      CALL timestop(0.0E0_dp,handle)

 END SUBROUTINE mc_moves_release

! *****************************************************************************
!!****s* mc_move_control.F/q_move_reinit
!!
!!   NAME
!!     q_move_reinit
!!
!!   SYNOPSIS
!!     Subroutine q_move_reinit(moves)
!!           Type(mc_moves_type), Pointer          :: moves
!!     End Subroutine q_move_reinit
!!
!!   FUNCTION
!!     sets all qsuccess counters back to zero
!!
!!   ARGUMENTS
!!     - moves: the move structure to update
!!  
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE move_q_reinit ( moves )

    TYPE(mc_moves_type), POINTER             :: moves

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "move_q_reinit"

    INTEGER                                  :: handle

!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

    CALL timeset(routineN,'I','',handle)

! set all the counters equal to zero
      moves%bond%qsuccesses=0
      moves%angle%qsuccesses=0
      moves%trans%qsuccesses=0
      moves%rot%qsuccesses=0
      moves%volume%qsuccesses=0
      moves%qtrans_dis=0.0E0_dp

 ! end the timing
      CALL timestop(0.0E0_dp,handle)

 END SUBROUTINE move_q_reinit

! *****************************************************************************
!!****s* mc_move_control.F/q_move_accept
!!
!!   NAME
!!     q_move_accept
!!
!!   FUNCTION
!!     updates accepted moves in the given structure...assumes you've been
!!     recording all successful moves in "qsuccesses"...this was done to
!!     compensate for doing multiple inner moves between Quickstep moves 
!!     (which determine ultimate acceptance of moves)
!!
!!   ARGUMENTS
!!     - moves: the move structure to update
!!  
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE q_move_accept(moves)

    TYPE(mc_moves_type), POINTER             :: moves

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "q_move_accept"

    INTEGER                                  :: handle

!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

    CALL timeset(routineN,'I','',handle)

! change the number of successful moves for the total move counter
    moves%bond%successes=moves%bond%successes&
           +moves%bond%qsuccesses
    moves%angle%successes=moves%angle%successes&
           +moves%angle%qsuccesses
    moves%trans%successes=moves%trans%successes&
           +moves%trans%qsuccesses
    moves%trans_dis=moves%trans_dis+moves%qtrans_dis
    moves%rot%successes=moves%rot%successes&
           +moves%rot%qsuccesses
    moves%volume%successes=moves%volume%successes&
           +moves%volume%qsuccesses

 ! change the number of successful moves for the move update counter
!    move_updates%bond%successes=move_updates%bond%successes&
!           +move_updates%bond%qsuccesses
!    move_updates%angle%successes=move_updates%angle%successes&
!           +move_updates%angle%qsuccesses
!    move_updates%trans%successes=move_updates%trans%successes&
!           +move_updates%trans%qsuccesses
!    move_updates%rot%successes=move_updates%rot%successes&
!           +move_updates%rot%qsuccesses
!    move_updates%volume%successes=move_updates%volume%successes&
!           +move_updates%volume%qsuccesses

! end the timing
      CALL timestop(0.0E0_dp,handle)

 END SUBROUTINE q_move_accept

! *****************************************************************************
!!****s* mc_move_control.F/q_move_accept
!!
!!   NAME
!!     q_move_accept
!!
!!   FUNCTION
!!     updates accepted moves in the given structure...assumes you've been
!!     recording all successful moves in "qsuccesses"...this was done to
!!     compensate for doing multiple inner moves between Quickstep moves 
!!     (which determine ultimate acceptance of moves)
!!
!!   ARGUMENTS
!!     - moves: the move structure to update
!!  
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE write_move_stats(moves,move_updates,nnstep,unit)

    TYPE(mc_moves_type), POINTER             :: moves, move_updates
    INTEGER, INTENT(IN)                      :: nnstep, unit

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "write_move_stats"

    INTEGER                                  :: handle

!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

    CALL timeset(routineN,'I','',handle)

    WRITE(unit,'(I6,2X,A,2X,I6,2X,I6)') nnstep,' classical_bond      ',&
      move_updates%classical_bond%successes+&
      move_updates%bond%successes,&
      move_updates%classical_bond%attempts+move_updates%bond%attempts
    WRITE(unit,'(I6,2X,A,2X,I6,2X,I6)') nnstep,' classical_angle     ',&
      move_updates%classical_angle%successes+&
      move_updates%angle%successes,&
      move_updates%classical_angle%attempts+move_updates%angle%attempts
    WRITE(unit,'(I6,2X,A,2X,I6,2X,I6)') nnstep,' classical_trans     ',&
      move_updates%classical_trans%successes+&
      move_updates%trans%successes,&
      move_updates%classical_trans%attempts+move_updates%trans%attempts
    WRITE(unit,'(I6,2X,A,2X,I6,2X,I6)') nnstep,' classical_rot       ',&
      move_updates%classical_rot%successes+move_updates%rot%successes,&
      move_updates%classical_rot%attempts+move_updates%rot%attempts
    WRITE(unit,'(I6,2X,A,2X,I6,2X,I6)') nnstep,' classical_volume    ',&
      moves%classical_volume%successes,&
      moves%classical_volume%attempts
    WRITE(unit,'(I6,2X,A,2X,I6,2X,I6)') nnstep,' bond      ',&
      moves%bond%successes,moves%bond%attempts
    WRITE(unit,'(I6,2X,A,2X,I6,2X,I6)') nnstep,' angle     ',&
      moves%angle%successes,moves%angle%attempts
    WRITE(unit,'(I6,2X,A,2X,I6,2X,I6)') nnstep,' trans     ',&
      moves%trans%successes,moves%trans%attempts
    WRITE(unit,'(I6,2X,A,2X,I6,2X,I6)') nnstep,' rot       ',&
      moves%rot%successes,moves%rot%attempts
    WRITE(unit,'(I6,2X,A,2X,I6,2X,I6)') nnstep,' swap      ',&
      moves%swap%successes,moves%swap%attempts
    WRITE(unit,'(I6,2X,A,2X,I6)') nnstep,' grown     ',&
      moves%grown
    WRITE(unit,'(I6,2X,A,2X,I6)') nnstep,' empty     ',&
      moves%empty
    WRITE(unit,'(I6,2X,A,2X,I6,2X,I6)') nnstep,' volume    ',&
      moves%volume%successes,moves%volume%attempts
    WRITE(unit,'(I6,2X,A,2X,I6,2X,I6)') nnstep,' Quickstep ',&
      moves%quickstep%successes,moves%quickstep%attempts

! end the timing
      CALL timestop(0.0E0_dp,handle)

 END SUBROUTINE write_move_stats

! *****************************************************************************
!!****s* mc_move_control.F/mc_move_update
!!
!!   NAME
!!     mc_move_update
!!
!!   FUNCTION
!!     updates the maximum displacements of a Monte Carlo simulation,
!!     based on the ratio of successful moves to attempts...tries to hit a
!!     target of 0.5 acceptance ratio
!!
!!   ARGUMENTS
!!     - mc_par: the mc parameters for the force env
!!     - move_updates: holds the accepted/attempted moves since the last
!!             update (or start of simulation)
!!     - nnstep: how many steps the simulation has run
!!     - flag: indicates which displacements to update..."volume" is for
!!             volume moves and "trans" is for everything else
!!  
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
 SUBROUTINE mc_move_update ( mc_par,move_updates,flag,nnstep )

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(mc_moves_type), POINTER             :: move_updates
    CHARACTER(LEN=*), INTENT(IN)             :: flag
    INTEGER, INTENT(IN)                      :: nnstep

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "mc_move_update"

    INTEGER                                  :: handle,rm
    REAL(KIND=dp)                            :: test_ratio,rmbond,rmangle,&
      rmtrans,rmrot,rmvolume
     
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! grab some stuff from mc_par
      CALL get_mc_par(mc_par,rmbond=rmbond,rmangle=rmangle,rmrot=rmrot,&
         rmtrans=rmtrans,rmvolume=rmvolume,rm=rm)

      SELECT CASE (flag)
      CASE DEFAULT
         WRITE(6,*) 'flag =',flag
         CALL stop_program( "mc_move_update","Wrong option passed")
      CASE ("trans")

! update the maximum displacement for bond length change
         IF( move_updates%bond%attempts .GT. 0 ) THEN

! first account for the extreme cases
            IF ( move_updates%bond%successes == 0 ) THEN
               rmbond=rmbond/2.0E0_dp
            ELSEIF ( move_updates%bond%successes == &
                     move_updates%bond%attempts ) THEN
                    rmbond=rmbond*2.0E0_dp
            ELSE
! now for the middle case
               test_ratio=REAL(move_updates%bond%successes,dp)&
                  /REAL(move_updates%bond%attempts,dp)/0.5E0_dp
               IF (test_ratio .GT. 2.0E0_dp) test_ratio=2.0E0_dp
               IF (test_ratio .LT. 0.5E0_dp) test_ratio=0.5E0_dp
               rmbond=rmbond*test_ratio
            ENDIF

! update and clear the counters
            move_updates%classical_bond%attempts=&
                move_updates%classical_bond%attempts+&
                move_updates%bond%attempts
            move_updates%classical_bond%successes=&
                move_updates%classical_bond%successes+&
                move_updates%bond%successes
            move_updates%bond%attempts=0
            move_updates%bond%successes=0

! write the new displacement to a file
            WRITE(rm,*) nnstep,' rmbond = ',&
                    rmbond*angstrom,' angstroms'

         ENDIF

! update the maximum displacement for bond length change
         IF( move_updates%angle%attempts .GT. 0 ) THEN

! first account for the extreme cases
            IF ( move_updates%angle%successes == 0 ) THEN
               rmangle=rmangle/2.0E0_dp
            ELSEIF ( move_updates%angle%successes == &
               move_updates%angle%attempts ) THEN
               rmangle=rmangle*2.0E0_dp
            ELSE
! now for the middle case
               test_ratio=REAL(move_updates%angle%successes,dp)&
                  /REAL(move_updates%angle%attempts,dp)/0.5E0_dp
               IF (test_ratio .GT. 2.0E0_dp) test_ratio=2.0E0_dp
               IF (test_ratio .LT. 0.5E0_dp) test_ratio=0.5E0_dp
               rmangle=rmangle*test_ratio
            ENDIF

! more than pi rotation is meaningless
               IF (rmangle .GT. pi) rmangle=pi

! clear the counters
            move_updates%classical_angle%attempts=&
                move_updates%classical_angle%attempts+&
                move_updates%angle%attempts
            move_updates%classical_angle%successes=&
                move_updates%classical_angle%successes+&
                move_updates%angle%successes
            move_updates%angle%attempts=0
            move_updates%angle%successes=0

! write the new displacement to a file
            WRITE(rm,*) nnstep,' rmangle = ',&
                    rmangle/pi*180.0E0_dp,' degrees'
         ENDIF

! update the maximum displacement for molecule translation
         IF( move_updates%trans%attempts .GT. 0 ) THEN

! first account for the extreme cases
            IF ( move_updates%trans%successes == 0 ) THEN
               rmtrans=rmtrans/2.0E0_dp
            ELSEIF ( move_updates%trans%successes == &
               move_updates%trans%attempts ) THEN
               rmtrans=rmtrans*2.0E0_dp
            ELSE
! now for the middle case
               test_ratio=REAL(move_updates%trans%successes,dp)&
                  /REAL(move_updates%trans%attempts,dp)/0.5E0_dp
               IF (test_ratio .GT. 2.0E0_dp) test_ratio=2.0E0_dp
               IF (test_ratio .LT. 0.5E0_dp) test_ratio=0.5E0_dp
               rmtrans=rmtrans*test_ratio
            ENDIF

 ! make an upper bound...10 a.u.
           IF (rmtrans .GT. 10.0E0_dp) rmtrans = 10.0E0_dp

 ! clear the counters
            move_updates%classical_trans%attempts=&
                move_updates%classical_trans%attempts+&
                move_updates%trans%attempts
            move_updates%classical_trans%successes=&
                move_updates%classical_trans%successes+&
                move_updates%trans%successes
            move_updates%trans%attempts=0
            move_updates%trans%successes=0

! write the new displacement to a file
            WRITE(rm,*) nnstep,' rmtrans = ',&
                    rmtrans*angstrom,' angstroms'
         ENDIF

! update the maximum displacement for molecule rotation
         IF( move_updates%rot%attempts .GT. 0 ) THEN

! first account for the extreme cases
            IF ( move_updates%rot%successes == 0 ) THEN
               rmrot=rmrot/2.0E0_dp

               IF (rmrot .GT. pi) rmrot=pi

            ELSEIF ( move_updates%rot%successes == &
                     move_updates%rot%attempts ) THEN
                 rmrot=rmrot/2.0E0_dp

! more than pi rotation is meaningless
               IF (rmrot .GT. pi) rmrot=pi

            ELSE
! now for the middle case
               test_ratio=REAL(move_updates%rot%successes,dp)&
                  /REAL(move_updates%rot%attempts,dp)/0.5E0_dp
               IF (test_ratio .GT. 2.0E0_dp) test_ratio=2.0E0_dp
               IF (test_ratio .LT. 0.5E0_dp) test_ratio=0.5E0_dp
               rmrot=rmrot*test_ratio

! more than pi rotation is meaningless
               IF (rmrot .GT. pi) rmrot=pi

            ENDIF

! clear the counters
            move_updates%classical_rot%attempts=&
                move_updates%classical_rot%attempts+&
                move_updates%rot%attempts
            move_updates%classical_rot%successes=&
                move_updates%classical_rot%successes+&
                move_updates%rot%successes
            move_updates%rot%attempts=0
            move_updates%rot%successes=0

! write the new displacement to a file
            WRITE(rm,*) nnstep,' rmrot = ',&
                    rmrot/pi*180.0E0_dp,' degrees'
         ENDIF
            
      CASE ("volume")

! update the maximum displacement for volume displacement
         IF ( move_updates%classical_volume%attempts .NE. 0) THEN

! first account for the extreme cases
            IF ( move_updates%classical_volume%successes == 0 ) THEN
               rmvolume=rmvolume/2.0E0_dp
               
            ELSEIF ( move_updates%classical_volume%successes == &
                  move_updates%classical_volume%attempts ) THEN
               rmvolume=rmvolume*2.0E0_dp
            ELSE
! now for the middle case
               test_ratio=REAL(move_updates%classical_volume%successes,dp)/&
                  REAL(move_updates%classical_volume%attempts,dp)/0.5E0_dp
               IF (test_ratio .GT. 2.0E0_dp) test_ratio=2.0E0_dp
               IF (test_ratio .LT. 0.5E0_dp) test_ratio=0.5E0_dp
               rmvolume=rmvolume*test_ratio

            ENDIF

! clear the counters
            move_updates%classical_volume%attempts=0
            move_updates%classical_volume%successes=0

! write the new displacement to a file
            WRITE(rm,*) nnstep,' rmvolume = ',&
                    rmvolume*angstrom**3,' angstroms^3'

         ENDIF
      END SELECT

! end the timing
      CALL timestop(0.0E0_dp,handle)

! grab some stuff from mc_par
      CALL set_mc_par(mc_par,rmbond=rmbond,rmangle=rmangle,rmrot=rmrot,&
         rmtrans=rmtrans,rmvolume=rmvolume)


  END SUBROUTINE mc_move_update

! -----------------------------------------------------------------------------

!**********************************************************************

END MODULE mc_move_control

!**********************************************************************

