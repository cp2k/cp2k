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
  USE mc_control,                      ONLY: mc_simulation_parameters_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_move_control"

  PUBLIC :: mc_moves_type, init_mc_moves, accattempt,&
            mc_move_update
 
  TYPE mc_moves_type
     TYPE ( accattempt ), POINTER :: transion
     TYPE ( accattempt ), POINTER :: trans
     TYPE ( accattempt ), POINTER :: rot
     TYPE ( accattempt ), POINTER :: volume
  END TYPE mc_moves_type

  TYPE accattempt
      INTEGER :: successes
      INTEGER :: attempts
  END TYPE accattempt

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE init_mc_moves ( moves )

!   Purpose: Canonical ensemble (NVT) Monte Carlo run using QUICKSTEP.

!   History: - Creation (10.14.2003)

!   ***************************************************************************

    TYPE(mc_moves_type), POINTER :: moves
    TYPE ( accattempt ),POINTER :: transion
    TYPE ( accattempt ),POINTER :: trans
    TYPE ( accattempt ),POINTER :: rot
    TYPE ( accattempt ),POINTER :: volume

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routineN = "init_mc_moves"

!   *** Local variables ***


!   *** Local Counters ***
      INTEGER :: stat,handle

!   ---------------------------------------------------------------------------

! begin the timing of the subroutine
    CALL timeset(routineN,'I','',handle)

! allocate all the structures...not sure why, but it won't work otherwise
      ALLOCATE (moves,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves",0)
      ALLOCATE (moves%transion,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%transion",0)
      ALLOCATE (moves%trans,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%trans",0)
      ALLOCATE (moves%rot,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%rot",0)
      ALLOCATE (moves%volume,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%volume",0)

! set all the counters equal to zero
      moves%transion%attempts=0
      moves%transion%successes=0
      moves%trans%attempts=0
      moves%trans%successes=0
      moves%rot%attempts=0
      moves%rot%successes=0
      moves%volume%attempts=0
      moves%volume%successes=0

 ! end the timing
      CALL timestop(0.0d0,handle)

 END SUBROUTINE init_mc_moves

 SUBROUTINE mc_move_update ( mc_par,moves,flag )

!   Purpose: Canonical ensemble (NVT) Monte Carlo run using QUICKSTEP.

!   History: - Creation (10.14.2003)

!   ***************************************************************************

    TYPE(mc_moves_type), POINTER :: moves
    TYPE ( mc_simulation_parameters_type ),POINTER  :: mc_par
    CHARACTER(LEN=*), INTENT(IN) :: flag

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_move_update"

!   *** Local variables ***


!   *** Local Counters ***
    INTEGER :: handle

!   ---------------------------------------------------------------------------

! begin the timing of the subroutine
      CALL timeset(routineN,'I','',handle)

      SELECT CASE (flag)
      CASE DEFAULT
         WRITE(6,*) 'flag =',flag
         CALL stop_program( "mc_move_update","Wrong option passed")
      CASE ("transion")

! update the maximum displacement for ion translation

! first account for the extreme cases
         IF ( moves%transion%successes == 0 ) THEN
            mc_par%rmion=mc_par%rmion*2.0d0
         ELSEIF ( moves%transion%successes == &
                 moves%transion%attempts ) THEN
            mc_par%rmion=mc_par%rmion/2.0d0
         ELSE
! now for the middle case
            mc_par%rmion=mc_par%rmion*((real(moves%transion%successes)&
                  /real(moves%transion%attempts))/0.5d0)
         ENDIF

! clear the counters
         moves%transion%attempts=0
         moves%transion%successes=0

      CASE ("trans")

! update the maximum displacement for ion translation

! first account for the extreme cases
         IF ( moves%trans%successes == 0 ) THEN
            mc_par%rmtrans=mc_par%rmtrans*2.0d0
         ELSEIF ( moves%trans%successes == &
                 moves%trans%attempts ) THEN
            mc_par%rmtrans=mc_par%rmtrans/2.0d0
         ELSE
! now for the middle case
            mc_par%rmtrans=mc_par%rmtrans*((real(moves%trans%successes)&
                  /real(moves%trans%attempts))/0.5d0)
         ENDIF

 ! clear the counters
         moves%trans%attempts=0
         moves%trans%successes=0
     

      CASE ("rot")

! update the maximum displacement for ion translation

! first account for the extreme cases
         IF ( moves%rot%successes == 0 ) THEN
            mc_par%rmrot=mc_par%rmrot*2.0d0

! more than pi rotation is meaningless
            IF (mc_par%rmrot .GT. 3.14159265358979d0) &
                 mc_par%rmrot=3.14159265358979d0

         ELSEIF ( moves%rot%successes == &
                 moves%rot%attempts ) THEN
            mc_par%rmrot=mc_par%rmrot/2.0d0
         ELSE
! now for the middle case
            mc_par%rmrot=mc_par%rmrot*((real(moves%rot%successes)&
                  /real(moves%rot%attempts))/0.5d0)

! more than pi rotation is meaningless
            IF (mc_par%rmrot .GT. 3.14159265358979d0) &
                 mc_par%rmrot=3.14159265358979d0

         ENDIF

! clear the counters
         moves%rot%attempts=0
         moves%rot%successes=0
         
     CASE ("volume")

! update the maximum displacement for volume displacement

! first account for the extreme cases
         IF ( moves%volume%successes == 0 ) THEN
            mc_par%rmvolume=mc_par%rmvolume*2.0d0

         ELSEIF ( moves%volume%successes == &
                 moves%volume%attempts ) THEN
            mc_par%rmvolume=mc_par%rmvolume/2.0d0
         ELSE
! now for the middle case
            mc_par%rmvolume=mc_par%rmvolume*((&
                   real(moves%volume%successes)&
                  /real(moves%volume%attempts))/0.5d0)
         ENDIF

! clear the counters
         moves%volume%attempts=0
         moves%volume%successes=0

      END SELECT

! end the timing
      CALL timestop(0.0d0,handle)

  END SUBROUTINE mc_move_update

! -----------------------------------------------------------------------------

!**********************************************************************

END MODULE mc_move_control

!**********************************************************************
