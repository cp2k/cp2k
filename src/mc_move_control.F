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
  USE force_env_types,                 ONLY: force_env_calculate_energy,&
                                             force_env_get,&
                                             force_env_type
  USE kinds,                           ONLY: dbl,&
                                             dp
  USE mc_control,                      ONLY: mc_simulation_parameters_type
  USE mathconstants,                   ONLY: pi
  USE physcon,                         ONLY: angstrom
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_move_control"

  PUBLIC :: mc_moves_type, mc_moves_p_type,init_mc_moves, accattempt,&
           mc_move_update,move_q_reinit,q_move_accept,mc_moves_release,&
           write_move_stats
 
! *** Types ***
  TYPE mc_moves_type
      TYPE ( accattempt ), POINTER :: classical_bond
      TYPE ( accattempt ), POINTER :: classical_angle
      TYPE ( accattempt ), POINTER :: classical_trans
      TYPE ( accattempt ), POINTER :: classical_rot
      TYPE ( accattempt ), POINTER :: bond
      TYPE ( accattempt ), POINTER :: angle
      TYPE ( accattempt ), POINTER :: trans
      TYPE ( accattempt ), POINTER :: rot
      TYPE ( accattempt ), POINTER :: swap
      TYPE ( accattempt ), POINTER :: volume
      TYPE ( accattempt ), POINTER :: Quickstep
      REAL(KIND = dp) :: classical_energy,trans_dis,qtrans_dis
      integer :: empty,grown
  END TYPE mc_moves_type

  TYPE accattempt
      INTEGER :: successes
      INTEGER :: qsuccesses
      INTEGER :: attempts
  END TYPE accattempt

  TYPE mc_moves_p_type
      TYPE(mc_moves_type), POINTER :: moves
  END TYPE mc_moves_p_type

! *****************************************************************************

CONTAINS

! *****************************************************************************
!!****s* mc_move_control.F/init_mc_moves
!!
!!   NAME
!!     init_mc_moves
!!
!!   SYNOPSIS
!!     Subroutine init_mc_moves (moves)
!!           Type(mc_moves_type), Pointer          :: moves
!!     End Subroutine init_mc_moves
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

    CHARACTER(LEN=*), PARAMETER    :: routineN = "init_mc_moves"

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
      moves%classical_energy=0.0d0
      moves%trans_dis=0.0d0
      moves%qtrans_dis=0.0d0
      moves%empty=0
      moves%grown=0
!      moves%force_create=1

 ! end the timing
      CALL timestop(0.0d0,handle)

 END SUBROUTINE init_mc_moves

! *****************************************************************************
!!****s* mc_move_control.F/mc_moves_release
!!
!!   NAME
!!     mc_moves_release
!!
!!   SYNOPSIS
!!     Subroutine mc_moves_release (moves)
!!           Type(mc_moves_type), Pointer          :: moves
!!     End Subroutine mc_moves_release
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

    CHARACTER(LEN=*), PARAMETER    :: routineN = "init_mc_moves"

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
      CALL timestop(0.0d0,handle)

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

    CHARACTER(LEN=*), PARAMETER     :: routineN = "init_mc_moves"

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
      moves%qtrans_dis=0.0d0

 ! end the timing
      CALL timestop(0.0d0,handle)

 END SUBROUTINE move_q_reinit

! *****************************************************************************
!!****s* mc_move_control.F/q_move_accept
!!
!!   NAME
!!     q_move_accept
!!
!!   SYNOPSIS
!!     Subroutine q_move_accept(moves)
!!           Type(mc_moves_type), Pointer          :: moves
!!     End Subroutine q_move_accept
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

    CHARACTER(LEN=*), PARAMETER       :: routineN = "init_mc_moves"

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
      CALL timestop(0.0d0,handle)

 END SUBROUTINE q_move_accept

! *****************************************************************************
!!****s* mc_move_control.F/q_move_accept
!!
!!   NAME
!!     q_move_accept
!!
!!   SYNOPSIS
!!     Subroutine q_move_accept(moves)
!!           Type(mc_moves_type), Pointer          :: moves
!!     End Subroutine q_move_accept
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

    TYPE(mc_moves_type), POINTER             :: moves,move_updates
    INTEGER,INTENT(IN)                       :: nnstep,unit

    CHARACTER(LEN=*), PARAMETER       :: routineN = "init_mc_moves"

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
      CALL timestop(0.0d0,handle)

 END SUBROUTINE write_move_stats

! *****************************************************************************
!!****s* mc_move_control.F/mc_move_update
!!
!!   NAME
!!     mc_move_update
!!
!!   SYNOPSIS
!!     Subroutine mc_move_update ( mc_par,move_updates,flag,nnstep )
!!           Type(mc_simulation_parameters_type),Pointer
!!                                                 :: mc_par
!!           Type(mc_moves_type), Pointer          :: move_updates
!!           Integer,Intent(IN)                    :: nnstep
!!           Character(Len=*),Intent(IN)           :: flag
!!     End Subroutine mc_move_update
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
 
    CHARACTER(LEN=*), PARAMETER       :: routineN = "mc_move_update"

    INTEGER                                  :: handle

    REAL(KIND=dp) :: test_ratio

!   ---------------------------------------------------------------------------

! begin the timing of the subroutine
      CALL timeset(routineN,'I','',handle)

      SELECT CASE (flag)
      CASE DEFAULT
         WRITE(6,*) 'flag =',flag
         CALL stop_program( "mc_move_update","Wrong option passed")
      CASE ("trans")

! update the maximum displacement for bond length change
         IF( move_updates%bond%attempts .GT. 0 ) THEN

! first account for the extreme cases
            IF ( move_updates%bond%successes == 0 ) THEN
               mc_par%rmbond=mc_par%rmbond/2.0d0
            ELSEIF ( move_updates%bond%successes == &
                     move_updates%bond%attempts ) THEN
                    mc_par%rmbond=mc_par%rmbond*2.0d0
            ELSE
! now for the middle case
               test_ratio=REAL(move_updates%bond%successes,dbl)&
                  /REAL(move_updates%bond%attempts,dbl)/0.5d0
               IF (test_ratio .GT. 2.0d0) test_ratio=2.0D0
               IF (test_ratio .LT. 0.5d0) test_ratio=0.5D0
               mc_par%rmbond=mc_par%rmbond*test_ratio
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
            WRITE(mc_par%rm,*) nnstep,' rmbond = ',&
                    mc_par%rmbond*angstrom,' angstroms'

         ENDIF

! update the maximum displacement for bond length change
         IF( move_updates%angle%attempts .GT. 0 ) THEN

! first account for the extreme cases
            IF ( move_updates%angle%successes == 0 ) THEN
               mc_par%rmangle=mc_par%rmangle/2.0d0
            ELSEIF ( move_updates%angle%successes == &
               move_updates%angle%attempts ) THEN
               mc_par%rmangle=mc_par%rmangle*2.0d0
            ELSE
! now for the middle case
               test_ratio=REAL(move_updates%angle%successes,dbl)&
                  /REAL(move_updates%angle%attempts,dbl)/0.5d0
               IF (test_ratio .GT. 2.0d0) test_ratio=2.0D0
               IF (test_ratio .LT. 0.5d0) test_ratio=0.5D0
               mc_par%rmangle=mc_par%rmangle*test_ratio
            ENDIF

! more than pi rotation is meaningless
               IF (mc_par%rmangle .GT. pi) &
                    mc_par%rmangle=pi

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
            WRITE(mc_par%rm,*) nnstep,' rmangle = ',&
                    mc_par%rmangle/pi*180.0d0,' degrees'
         ENDIF

! update the maximum displacement for molecule translation
         IF( move_updates%trans%attempts .GT. 0 ) THEN

! first account for the extreme cases
            IF ( move_updates%trans%successes == 0 ) THEN
               mc_par%rmtrans=mc_par%rmtrans/2.0d0
            ELSEIF ( move_updates%trans%successes == &
               move_updates%trans%attempts ) THEN
               mc_par%rmtrans=mc_par%rmtrans*2.0d0
            ELSE
! now for the middle case
               test_ratio=REAL(move_updates%trans%successes,dbl)&
                  /REAL(move_updates%trans%attempts,dbl)/0.5d0
               IF (test_ratio .GT. 2.0d0) test_ratio=2.0D0
               IF (test_ratio .LT. 0.5d0) test_ratio=0.5D0
               mc_par%rmtrans=mc_par%rmtrans*test_ratio
            ENDIF

 ! make an upper bound...10 a.u.
           IF (mc_par%rmtrans .GT. 10.0d0) mc_par%rmtrans = 10.0d0

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
            WRITE(mc_par%rm,*) nnstep,' rmtrans = ',&
                    mc_par%rmtrans*angstrom,' angstroms'
         ENDIF

! update the maximum displacement for molecule rotation
         IF( move_updates%rot%attempts .GT. 0 ) THEN

! first account for the extreme cases
            IF ( move_updates%rot%successes == 0 ) THEN
               mc_par%rmrot=mc_par%rmrot/2.0d0

               IF (mc_par%rmrot .GT. pi) &
                                  mc_par%rmrot=pi

            ELSEIF ( move_updates%rot%successes == &
                     move_updates%rot%attempts ) THEN
                 mc_par%rmrot=mc_par%rmrot/2.0d0

! more than pi rotation is meaningless
               IF (mc_par%rmrot .GT. pi) &
                    mc_par%rmrot=pi

            ELSE
! now for the middle case
               test_ratio=REAL(move_updates%rot%successes,dbl)&
                  /REAL(move_updates%rot%attempts,dbl)/0.5d0
               IF (test_ratio .GT. 2.0d0) test_ratio=2.0D0
               IF (test_ratio .LT. 0.5d0) test_ratio=0.5D0
               mc_par%rmrot=mc_par%rmrot*test_ratio

! more than pi rotation is meaningless
               IF (mc_par%rmrot .GT. pi) &
                    mc_par%rmrot=pi

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
            WRITE(mc_par%rm,*) nnstep,' rmrot = ',&
                    mc_par%rmrot/pi*180.0d0,' degrees'
         ENDIF
            
     CASE ("volume")

! update the maximum displacement for volume displacement

! first account for the extreme cases
         IF ( move_updates%volume%successes == 0 ) THEN
            mc_par%rmvolume=mc_par%rmvolume/2.0d0

         ELSEIF ( move_updates%volume%successes == &
                 move_updates%volume%attempts ) THEN
            mc_par%rmvolume=mc_par%rmvolume*2.0d0
         ELSE
! now for the middle case
            mc_par%rmvolume=mc_par%rmvolume*((&
                   REAL(move_updates%volume%successes,dbl)&
                  /REAL(move_updates%volume%attempts,dbl))/0.5d0)
         ENDIF

! clear the counters
         move_updates%volume%attempts=0
         move_updates%volume%successes=0

      END SELECT

! end the timing
      CALL timestop(0.0d0,handle)

  END SUBROUTINE mc_move_update

! -----------------------------------------------------------------------------

!**********************************************************************

END MODULE mc_move_control

!**********************************************************************
