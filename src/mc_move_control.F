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
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_move_control"

  PUBLIC :: mc_moves_type, init_mc_moves, accattempt,&
            mc_move_update, move_q_reinit, q_move_accept,&
            mc_energy_restart
 
  TYPE mc_moves_type
     TYPE ( accattempt ), POINTER :: bond
     TYPE ( accattempt ), POINTER :: angle
     TYPE ( accattempt ), POINTER :: trans
     TYPE ( accattempt ), POINTER :: rot
     TYPE ( accattempt ), POINTER :: volume
     TYPE ( accattempt ), POINTER :: Quickstep
     REAL(KIND = dp) :: classical_energy,trans_dis,qtrans_dis
!     INTEGER :: force_create
  END TYPE mc_moves_type

  TYPE accattempt
      INTEGER :: successes
      INTEGER :: qsuccesses
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

    TYPE(mc_moves_type), POINTER             :: moves

    CHARACTER(LEN=*), PARAMETER    :: routineN = "init_mc_moves"

    INTEGER                                  :: handle, stat

!   *** Local Counters ***
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
      ALLOCATE (moves%volume,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%volume",0)
      ALLOCATE (moves%Quickstep,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "moves%Quickstep",0)
! set all the counters equal to zero
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
      moves%Quickstep%attempts=0
      moves%Quickstep%successes=0
      moves%Quickstep%qsuccesses=0
      moves%classical_energy=0.0d0
      moves%trans_dis=0.0d0
      moves%qtrans_dis=0.0d0
!      moves%force_create=1

 ! end the timing
      CALL timestop(0.0d0,handle)

 END SUBROUTINE init_mc_moves

SUBROUTINE move_q_reinit ( moves )

!   Purpose: Canonical ensemble (NVT) Monte Carlo run using QUICKSTEP.

!   History: - Creation (10.14.2003)

!   ***************************************************************************

    TYPE(mc_moves_type), POINTER             :: moves

    CHARACTER(LEN=*), PARAMETER     :: routineN = "init_mc_moves"

    INTEGER                                  :: handle

!   *** Local Counters ***
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

SUBROUTINE q_move_accept(moves,move_updates)

!   Purpose: Canonical ensemble (NVT) Monte Carlo run using QUICKSTEP.

!   History: - Creation (10.14.2003)

!   ***************************************************************************

    TYPE(mc_moves_type), POINTER             :: moves, move_updates

    CHARACTER(LEN=*), PARAMETER       :: routineN = "init_mc_moves"

    INTEGER                                  :: handle

!   *** Local Counters ***
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
    move_updates%bond%successes=move_updates%bond%successes&
           +move_updates%bond%qsuccesses
    move_updates%angle%successes=move_updates%angle%successes&
           +move_updates%angle%qsuccesses
    move_updates%trans%successes=move_updates%trans%successes&
           +move_updates%trans%qsuccesses
    move_updates%rot%successes=move_updates%rot%successes&
           +move_updates%rot%qsuccesses
    move_updates%volume%successes=move_updates%volume%successes&
           +move_updates%volume%qsuccesses

! end the timing
      CALL timestop(0.0d0,handle)

 END SUBROUTINE q_move_accept
!***************************************************************************

!***************************************************************************
 SUBROUTINE mc_move_update ( mc_par,move_updates,flag )

!   Purpose: Canonical ensemble (NVT) Monte Carlo run using QUICKSTEP.

!   History: - Creation (10.14.2003)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(mc_moves_type), POINTER             :: move_updates
    CHARACTER(LEN=*), INTENT(IN)             :: flag

    CHARACTER(LEN=*), PARAMETER       :: routineN = "mc_move_update"

    INTEGER                                  :: handle

!   *** Local Counters ***
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

      SELECT CASE (flag)
      CASE DEFAULT
         WRITE(6,*) 'flag =',flag
         CALL stop_program( "mc_move_update","Wrong option passed")
      CASE ("bond")

! update the maximum displacement for bond length change

! first account for the extreme cases
         IF ( move_updates%bond%successes == 0 ) THEN
            mc_par%rmbond=mc_par%rmbond/2.0d0
         ELSEIF ( move_updates%bond%successes == &
                 move_updates%bond%attempts ) THEN
            mc_par%rmbond=mc_par%rmbond*2.0d0
         ELSE
! now for the middle case
            mc_par%rmbond=mc_par%rmbond*(&
                  (REAL(move_updates%bond%successes,dbl)&
                  /REAL(move_updates%bond%attempts,dbl))/0.5d0)
         ENDIF

! clear the counters
         move_updates%bond%attempts=0
         move_updates%bond%successes=0

      CASE ("angle")

! update the maximum displacement for bond length change

! first account for the extreme cases
         IF ( move_updates%angle%successes == 0 ) THEN
            mc_par%rmangle=mc_par%rmangle/2.0d0
         ELSEIF ( move_updates%angle%successes == &
                 move_updates%angle%attempts ) THEN
            mc_par%rmangle=mc_par%rmangle*2.0d0
         ELSE
! now for the middle case
            mc_par%rmangle=mc_par%rmangle*(&
                  (REAL(move_updates%angle%successes,dbl)&
                  /REAL(move_updates%angle%attempts,dbl))/0.5d0)
         ENDIF

! clear the counters
         move_updates%angle%attempts=0
         move_updates%angle%successes=0

      CASE ("trans")

! update the maximum displacement for molecule translation

! first account for the extreme cases
         IF ( move_updates%trans%successes == 0 ) THEN
            mc_par%rmtrans=mc_par%rmtrans/2.0d0
         ELSEIF ( move_updates%trans%successes == &
                 move_updates%trans%attempts ) THEN
            mc_par%rmtrans=mc_par%rmtrans*2.0d0
         ELSE
! now for the middle case
            mc_par%rmtrans=mc_par%rmtrans*((&
                  REAL(move_updates%trans%successes,dbl)&
                  /REAL(move_updates%trans%attempts,dbl))/0.5d0)
         ENDIF

 ! clear the counters
         move_updates%trans%attempts=0
         move_updates%trans%successes=0
     

      CASE ("rot")

! update the maximum displacement for molecule rotation

! first account for the extreme cases
         IF ( move_updates%rot%successes == 0 ) THEN
            mc_par%rmrot=mc_par%rmrot/2.0d0

            IF (mc_par%rmrot .GT. 3.14159265358979d0) &
                 mc_par%rmrot=3.14159265358979d0

         ELSEIF ( move_updates%rot%successes == &
                 move_updates%rot%attempts ) THEN
            mc_par%rmrot=mc_par%rmrot/2.0d0

! more than pi rotation is meaningless
            IF (mc_par%rmrot .GT. 3.14159265358979d0) &
                 mc_par%rmrot=3.14159265358979d0

         ELSE
! now for the middle case
            mc_par%rmrot=mc_par%rmrot*((&
                  REAL(move_updates%rot%successes,dbl)&
                  /REAL(move_updates%rot%attempts,dbl))/0.5d0)

! more than pi rotation is meaningless
            IF (mc_par%rmrot .GT. 3.14159265358979d0) &
                 mc_par%rmrot=3.14159265358979d0

         ENDIF

! clear the counters
         move_updates%rot%attempts=0
         move_updates%rot%successes=0
         
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

 SUBROUTINE mc_energy_restart(force_env,energy,error)
! the system call might cause problems on various platforms, so
! we do it here
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), INTENT(OUT)               :: energy
    TYPE(cp_error_type), INTENT(INOUT), OPTIONAL   :: error


      CALL force_env_calculate_energy(force_env)
      CALL force_env_get(force_env,&
          potential_energy=energy,error=error)

  END SUBROUTINE mc_energy_restart

!**********************************************************************

END MODULE mc_move_control

!**********************************************************************
