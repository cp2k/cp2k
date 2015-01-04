!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief control the handling of the move data in Monte Carlo (MC) simulations
!> \par History
!>      none
!> \author Matthew J. McGrath  (10.16.2003)
! *****************************************************************************
MODULE mc_move_control
  
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE mc_types,                        ONLY: get_mc_molecule_info,&
                                             get_mc_par,&
                                             mc_molecule_info_type,&
                                             mc_moves_type,&
                                             mc_simpar_type,&
                                             set_mc_par
  USE physcon,                         ONLY: angstrom
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mc_move_control'

  PUBLIC :: init_mc_moves, &
           mc_move_update,move_q_reinit,q_move_accept,mc_moves_release,&
           write_move_stats

CONTAINS

! *****************************************************************************
!> \brief allocates and initializes the structure to record all move
!>      attempts/successes
!> \param moves the move structure to update
!>
!>    Suitable for parallel.
!> \author MJM
! *****************************************************************************
  SUBROUTINE init_mc_moves ( moves )

    TYPE(mc_moves_type), POINTER             :: moves

    CHARACTER(len=*), PARAMETER :: routineN = 'init_mc_moves', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat

! begin the timing of the subroutine

    CALL timeset(routineN,handle)

! allocate all the structures
      ALLOCATE (moves,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves",0)
      ALLOCATE (moves%bond,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%bond",0)
      ALLOCATE (moves%angle,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%angle",0)
      ALLOCATE (moves%dihedral,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%dihedral",0)
      ALLOCATE (moves%trans,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%trans",0)
      ALLOCATE (moves%rot,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%rot",0)
      ALLOCATE (moves%bias_bond,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%bias_bond",0)
      ALLOCATE (moves%bias_angle,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%bias_angle",0)
      ALLOCATE (moves%bias_dihedral,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%bias_dihedral",0)
      ALLOCATE (moves%bias_trans,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%bias_trans",0)
      ALLOCATE (moves%bias_rot,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%bias_rot",0)
      ALLOCATE (moves%volume,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%volume",0)
      ALLOCATE (moves%hmc,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%hmc",0)
      ALLOCATE (moves%avbmc_inin,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%avbmc_inin",0)
      ALLOCATE (moves%avbmc_inout,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%avbmc_inout",0)
      ALLOCATE (moves%avbmc_outin,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%avbmc_outin",0)
      ALLOCATE (moves%avbmc_outout,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%avbmc_outout",0)
      ALLOCATE (moves%swap,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%swap",0)
      ALLOCATE (moves%Quickstep,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%Quickstep",0)
! set all the counters equal to zero
      moves%bias_bond%attempts=0
      moves%bias_bond%successes=0
      moves%bias_bond%qsuccesses=0
      moves%bias_angle%attempts=0
      moves%bias_angle%successes=0
      moves%bias_angle%qsuccesses=0
      moves%bias_dihedral%attempts=0
      moves%bias_dihedral%successes=0
      moves%bias_dihedral%qsuccesses=0
      moves%bias_trans%attempts=0
      moves%bias_trans%successes=0
      moves%bias_trans%qsuccesses=0
      moves%bias_rot%attempts=0
      moves%bias_rot%successes=0
      moves%bias_rot%qsuccesses=0

      moves%bond%attempts=0
      moves%bond%successes=0
      moves%bond%qsuccesses=0
      moves%angle%attempts=0
      moves%angle%successes=0
      moves%angle%qsuccesses=0
      moves%dihedral%attempts=0
      moves%dihedral%successes=0
      moves%dihedral%qsuccesses=0
      moves%trans%attempts=0
      moves%trans%successes=0
      moves%trans%qsuccesses=0
      moves%rot%attempts=0
      moves%rot%successes=0
      moves%rot%qsuccesses=0
      moves%avbmc_inin%attempts=0
      moves%avbmc_inin%successes=0
      moves%avbmc_inin%qsuccesses=0
      moves%avbmc_inout%attempts=0
      moves%avbmc_inout%successes=0
      moves%avbmc_inout%qsuccesses=0
      moves%avbmc_outin%attempts=0
      moves%avbmc_outin%successes=0
      moves%avbmc_outin%qsuccesses=0
      moves%avbmc_outout%attempts=0
      moves%avbmc_outout%successes=0
      moves%avbmc_outout%qsuccesses=0
      moves%volume%attempts=0
      moves%volume%successes=0
      moves%volume%qsuccesses=0
      moves%hmc%attempts=0
      moves%hmc%successes=0
      moves%hmc%qsuccesses=0
      moves%swap%attempts=0
      moves%swap%successes=0
      moves%swap%qsuccesses=0
      moves%Quickstep%attempts=0
      moves%Quickstep%successes=0
      moves%Quickstep%qsuccesses=0
      moves%trans_dis=0.0E0_dp
      moves%qtrans_dis=0.0E0_dp
      moves%empty=0
      moves%empty_conf=0
      moves%empty_avbmc=0
      moves%grown=0
!      moves%force_create=1

 ! end the timing
      CALL timestop(handle)

 END SUBROUTINE init_mc_moves

! *****************************************************************************
!> \brief deallocates all the structures and nullifies the pointer
!> \param moves the move structure to release
!>
!>    Suitable for parallel.
!> \author MJM
! *****************************************************************************
  SUBROUTINE mc_moves_release ( moves )

    TYPE(mc_moves_type), POINTER             :: moves

    CHARACTER(len=*), PARAMETER :: routineN = 'mc_moves_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat

! begin the timing of the subroutine

    CALL timeset(routineN,handle)

! allocate all the structures
      DEALLOCATE (moves%bond,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%bond",0)
      DEALLOCATE (moves%angle,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%angle",0)
      DEALLOCATE (moves%dihedral,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%dihedral",0)
      DEALLOCATE (moves%trans,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%trans",0)
      DEALLOCATE (moves%rot,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%rot",0)
      DEALLOCATE (moves%bias_bond,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%bias_bond",0)
      DEALLOCATE (moves%bias_angle,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%bias_angle",0)
      DEALLOCATE (moves%bias_dihedral,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%bias_dihedral",0)
      DEALLOCATE (moves%bias_trans,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%bias_trans",0)
      DEALLOCATE (moves%bias_rot,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%bias_rot",0)
      DEALLOCATE (moves%volume,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%volume",0)
      DEALLOCATE (moves%hmc,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%hmc",0)
      DEALLOCATE (moves%avbmc_inin,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%avbmc_inin",0)
      DEALLOCATE (moves%avbmc_inout,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%avbmc_inout",0)
      DEALLOCATE (moves%avbmc_outin,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%avbmc_outin",0)
      DEALLOCATE (moves%avbmc_outout,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%avbmc_outout",0)
      DEALLOCATE (moves%swap,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%swap",0)
      DEALLOCATE (moves%Quickstep,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves%Quickstep",0)

      DEALLOCATE (moves,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "moves",0)

! now nullify the moves
      NULLIFY(moves)

 ! end the timing
      CALL timestop(handle)

 END SUBROUTINE mc_moves_release

! *****************************************************************************
!> \brief sets all qsuccess counters back to zero
!> \param moves the move structure to update
!> \param lbias are we biasing translations/rotations/conformational changes
!>        with a different potential?
!>
!>    Suitable for parallel.
!> \author MJM
! *****************************************************************************
SUBROUTINE move_q_reinit ( moves , lbias )

    TYPE(mc_moves_type), POINTER             :: moves
    LOGICAL, INTENT(IN)                      :: lbias

    CHARACTER(len=*), PARAMETER :: routineN = 'move_q_reinit', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

! begin the timing of the subroutine

      CALL timeset(routineN,handle)

! set all the counters equal to zero
      IF(lbias) THEN
         moves%bias_bond%qsuccesses=0
         moves%bias_angle%qsuccesses=0
         moves%bias_dihedral%qsuccesses=0
         moves%bias_trans%qsuccesses=0
         moves%bias_rot%qsuccesses=0
      ELSE
         moves%bond%qsuccesses=0
         moves%angle%qsuccesses=0
         moves%dihedral%qsuccesses=0
         moves%trans%qsuccesses=0
         moves%rot%qsuccesses=0
         moves%volume%qsuccesses=0
         moves%hmc%qsuccesses=0
         moves%qtrans_dis=0.0E0_dp
      ENDIF

 ! end the timing
      CALL timestop(handle)

 END SUBROUTINE move_q_reinit

! *****************************************************************************
!> \brief updates accepted moves in the given structure...assumes you've been
!>      recording all successful moves in "qsuccesses"...this was done to
!>      compensate for doing multiple inner moves between Quickstep moves
!>      (which determine ultimate acceptance of moves)
!> \param moves the move structure to update
!> \param lbias are we biasing non-swap particle moves with a cheaper potential
!>
!>    Suitable for parallel.
!> \author MJM
! *****************************************************************************
SUBROUTINE q_move_accept(moves,lbias)

    TYPE(mc_moves_type), POINTER             :: moves
    LOGICAL, INTENT(IN)                      :: lbias

    CHARACTER(len=*), PARAMETER :: routineN = 'q_move_accept', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

! begin the timing of the subroutine

      CALL timeset(routineN,handle)

      IF(lbias) THEN
! change the number of successful moves for the total move counter
         moves%bias_bond%successes=moves%bias_bond%successes&
            +moves%bias_bond%qsuccesses
         moves%bias_angle%successes=moves%bias_angle%successes&
            +moves%bias_angle%qsuccesses
         moves%bias_dihedral%successes=moves%bias_dihedral%successes&
            +moves%bias_dihedral%qsuccesses
         moves%bias_trans%successes=moves%bias_trans%successes&
            +moves%bias_trans%qsuccesses
         moves%bias_rot%successes=moves%bias_rot%successes&
            +moves%bias_rot%qsuccesses
      ELSE
! change the number of successful moves for the total move counter
         moves%bond%successes=moves%bond%successes&
            +moves%bond%qsuccesses
         moves%angle%successes=moves%angle%successes&
            +moves%angle%qsuccesses
         moves%dihedral%successes=moves%dihedral%successes&
            +moves%dihedral%qsuccesses
         moves%trans%successes=moves%trans%successes&
            +moves%trans%qsuccesses
         moves%rot%successes=moves%rot%successes&
            +moves%rot%qsuccesses
         moves%hmc%successes=moves%hmc%successes&
            +moves%hmc%qsuccesses
         moves%volume%successes=moves%volume%successes&
            +moves%volume%qsuccesses
         moves%avbmc_inin%successes=moves%avbmc_inin%successes&
            +moves%avbmc_inin%qsuccesses
         moves%avbmc_inout%successes=moves%avbmc_inout%successes&
            +moves%avbmc_inout%qsuccesses
         moves%avbmc_outin%successes=moves%avbmc_outin%successes&
            +moves%avbmc_outin%qsuccesses
         moves%avbmc_outout%successes=moves%avbmc_outout%successes&
            +moves%avbmc_outout%qsuccesses

         moves%trans_dis=moves%trans_dis+moves%qtrans_dis
      ENDIF

! end the timing
      CALL timestop(handle)

 END SUBROUTINE q_move_accept

! *****************************************************************************
!> \brief writes the number of accepted and attempted moves to a file for
!>      the various move types
!> \param moves the structure containing the move data
!> \param nnstep what step we're on
!> \param unit the unit of the file we're writing to
!>
!>    Use only in serial.
!> \author MJM
! *****************************************************************************
SUBROUTINE write_move_stats(moves,nnstep,unit)

    TYPE(mc_moves_type), POINTER             :: moves
    INTEGER, INTENT(IN)                      :: nnstep, unit

    CHARACTER(len=*), PARAMETER :: routineN = 'write_move_stats', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

! begin the timing of the subroutine

    CALL timeset(routineN,handle)

    WRITE(unit,1000) nnstep,' bias_bond      ',&
      moves%bias_bond%successes,moves%bias_bond%attempts
    WRITE(unit,1000) nnstep,' bias_angle      ',&
      moves%bias_angle%successes,moves%bias_angle%attempts
    WRITE(unit,1000) nnstep,' bias_dihedral      ',&
      moves%bias_dihedral%successes,moves%bias_dihedral%attempts
    WRITE(unit,1000) nnstep,' bias_trans      ',&
      moves%bias_trans%successes,moves%bias_trans%attempts
    WRITE(unit,1000) nnstep,' bias_rot      ',&
      moves%bias_rot%successes,moves%bias_rot%attempts

    WRITE(unit,1000) nnstep,' bond      ',&
      moves%bond%successes,moves%bond%attempts
    WRITE(unit,1000) nnstep,' angle     ',&
      moves%angle%successes,moves%angle%attempts
    WRITE(unit,1000) nnstep,' dihedral     ',&
      moves%dihedral%successes,moves%dihedral%attempts
    WRITE(unit,1000) nnstep,' trans     ',&
      moves%trans%successes,moves%trans%attempts
    WRITE(unit,1000) nnstep,' rot       ',&
      moves%rot%successes,moves%rot%attempts
    WRITE(unit,1000) nnstep,' swap      ',&
      moves%swap%successes,moves%swap%attempts
    WRITE(unit,1001) nnstep,' grown     ',&
      moves%grown
    WRITE(unit,1001) nnstep,' empty_swap     ',&
      moves%empty
    WRITE(unit,1001) nnstep,' empty_conf     ',&
      moves%empty_conf
    WRITE(unit,1000) nnstep,' volume    ',&
      moves%volume%successes,moves%volume%attempts
    WRITE(unit,1000) nnstep,' HMC    ',&
      moves%hmc%successes,moves%hmc%attempts
    WRITE(unit,1000) nnstep,' avbmc_inin  ',&
      moves%avbmc_inin%successes,moves%avbmc_inin%attempts
    WRITE(unit,1000) nnstep,' avbmc_inout  ',&
      moves%avbmc_inout%successes,moves%avbmc_inout%attempts
    WRITE(unit,1000) nnstep,' avbmc_outin  ',&
      moves%avbmc_outin%successes,moves%avbmc_outin%attempts
    WRITE(unit,1000) nnstep,' avbmc_outout ',&
      moves%avbmc_outout%successes,moves%avbmc_outout%attempts
    WRITE(unit,1001) nnstep,' empty_avbmc     ',&
      moves%empty_avbmc
    WRITE(unit,1000) nnstep,' Quickstep ',&
      moves%quickstep%successes,moves%quickstep%attempts

1000 FORMAT(I10,2X,A,2X,I10,2X,I10)
1001 FORMAT(I10,2X,A,2X,I10)
! end the timing
      CALL timestop(handle)

 END SUBROUTINE write_move_stats

! *****************************************************************************
!> \brief updates the maximum displacements of a Monte Carlo simulation,
!>      based on the ratio of successful moves to attempts...tries to hit a
!>      target of 0.5 acceptance ratio
!> \param mc_par the mc parameters for the force env
!> \param move_updates holds the accepted/attempted moves since the last
!>              update (or start of simulation)
!> \param molecule_type ...
!> \param flag indicates which displacements to update..."volume" is for
!>              volume moves and "trans" is for everything else
!> \param nnstep how many steps the simulation has run
!> \param ionode is this the main CPU running this job?
!>
!>    Suitable for parallel.
!> \author MJM
! *****************************************************************************
 SUBROUTINE mc_move_update ( mc_par,move_updates,molecule_type,flag,&
      nnstep,ionode )

    TYPE(mc_simpar_type), POINTER            :: mc_par
    TYPE(mc_moves_type), POINTER             :: move_updates
    INTEGER, INTENT(IN)                      :: molecule_type
    CHARACTER(LEN=*), INTENT(IN)             :: flag
    INTEGER, INTENT(IN)                      :: nnstep
    LOGICAL, INTENT(IN)                      :: ionode

    CHARACTER(len=*), PARAMETER :: routineN = 'mc_move_update', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nmol_types, rm
    REAL(dp), DIMENSION(:), POINTER          :: rmangle, rmbond, rmdihedral, &
                                                rmrot, rmtrans
    REAL(KIND=dp)                            :: rmvolume, test_ratio
    TYPE(mc_molecule_info_type), POINTER     :: mc_molecule_info

! begin the timing of the subroutine

      CALL timeset(routineN,handle)

      NULLIFY(rmangle,rmbond,rmdihedral,rmrot,rmtrans)

! grab some stuff from mc_par
      CALL get_mc_par(mc_par,rmbond=rmbond,rmangle=rmangle,rmrot=rmrot,&
         rmtrans=rmtrans,rmvolume=rmvolume,rm=rm,rmdihedral=rmdihedral,&
         mc_molecule_info=mc_molecule_info)
      CALL get_mc_molecule_info(mc_molecule_info,nmol_types=nmol_types)

      SELECT CASE (flag)
      CASE DEFAULT
         WRITE(*,*) 'flag =',flag
         CALL stop_program(routineN,moduleN,__LINE__,"Wrong option passed")
      CASE ("trans")

! we need to update all the displacements for every molecule type
         IF(ionode) WRITE(rm,*) nnstep,' Data for molecule type ',&
              molecule_type

! update the maximum displacement for bond length change
         IF( move_updates%bias_bond%attempts .GT. 0 ) THEN

! first account for the extreme cases
            IF ( move_updates%bias_bond%successes == 0 ) THEN
               rmbond(molecule_type)=rmbond(molecule_type)/2.0E0_dp
            ELSEIF ( move_updates%bias_bond%successes == &
                 move_updates%bias_bond%attempts ) THEN
               rmbond(molecule_type)=rmbond(molecule_type)*2.0E0_dp
            ELSE
! now for the middle case
               test_ratio=REAL(move_updates%bias_bond%successes,dp)&
                    /REAL(move_updates%bias_bond%attempts,dp)/0.5E0_dp
               IF (test_ratio .GT. 2.0E0_dp) test_ratio=2.0E0_dp
               IF (test_ratio .LT. 0.5E0_dp) test_ratio=0.5E0_dp
               rmbond(molecule_type)=rmbond(molecule_type)*test_ratio
            ENDIF

! update and clear the counters
            move_updates%bias_bond%attempts=0
            move_updates%bias_bond%successes=0

! write the new displacement to a file
            IF(ionode) WRITE(rm,*) nnstep,' rmbond = ',&
                 rmbond(molecule_type)*angstrom,' angstroms'

         ENDIF

! update the maximum displacement for bond angle change
         IF( move_updates%bias_angle%attempts .GT. 0 ) THEN

! first account for the extreme cases
            IF ( move_updates%bias_angle%successes == 0 ) THEN
               rmangle(molecule_type)=rmangle(molecule_type)/2.0E0_dp
            ELSEIF ( move_updates%bias_angle%successes == &
                 move_updates%bias_angle%attempts ) THEN
               rmangle(molecule_type)=rmangle(molecule_type)*2.0E0_dp
            ELSE
! now for the middle case
               test_ratio=REAL(move_updates%bias_angle%successes,dp)&
                    /REAL(move_updates%bias_angle%attempts,dp)/0.5E0_dp
               IF (test_ratio .GT. 2.0E0_dp) test_ratio=2.0E0_dp
               IF (test_ratio .LT. 0.5E0_dp) test_ratio=0.5E0_dp
               rmangle(molecule_type)=rmangle(molecule_type)*test_ratio
            ENDIF

! more than pi changes meaningless
            IF (rmangle(molecule_type) .GT. pi) rmangle(molecule_type)=pi

! clear the counters
            move_updates%bias_angle%attempts=0
            move_updates%bias_angle%successes=0

! write the new displacement to a file
            IF(ionode) WRITE(rm,*) nnstep,' rmangle = ',&
                 rmangle(molecule_type)/pi*180.0E0_dp,' degrees'
         ENDIF

! update the maximum displacement for a dihedral change
         IF( move_updates%bias_dihedral%attempts .GT. 0 ) THEN

! first account for the extreme cases
            IF ( move_updates%bias_dihedral%successes == 0 ) THEN
               rmdihedral(molecule_type)=rmdihedral(molecule_type)/2.0E0_dp
            ELSEIF ( move_updates%bias_dihedral%successes == &
                 move_updates%bias_dihedral%attempts ) THEN
               rmdihedral(molecule_type)=rmdihedral(molecule_type)*2.0E0_dp
            ELSE
! now for the middle case
               test_ratio=REAL(move_updates%bias_dihedral%successes,dp)&
                    /REAL(move_updates%bias_dihedral%attempts,dp)/0.5E0_dp
               IF (test_ratio .GT. 2.0E0_dp) test_ratio=2.0E0_dp
               IF (test_ratio .LT. 0.5E0_dp) test_ratio=0.5E0_dp
               rmdihedral(molecule_type)=rmdihedral(molecule_type)*test_ratio
            ENDIF

! more than pi changes meaningless
            IF (rmdihedral(molecule_type) .GT. pi) rmdihedral(molecule_type)=pi

! clear the counters
            move_updates%bias_dihedral%attempts=0
            move_updates%bias_dihedral%successes=0

! write the new displacement to a file
            IF(ionode) WRITE(rm,*) nnstep,' rmdihedral = ',&
                 rmdihedral(molecule_type)/pi*180.0E0_dp,' degrees'
         ENDIF

! update the maximum displacement for molecule translation
         IF( move_updates%bias_trans%attempts .GT. 0 ) THEN

! first account for the extreme cases
            IF ( move_updates%bias_trans%successes == 0 ) THEN
               rmtrans(molecule_type)=rmtrans(molecule_type)/2.0E0_dp
            ELSEIF ( move_updates%bias_trans%successes == &
                 move_updates%bias_trans%attempts ) THEN
               rmtrans(molecule_type)=rmtrans(molecule_type)*2.0E0_dp
            ELSE
! now for the middle case
               test_ratio=REAL(move_updates%bias_trans%successes,dp)&
                    /REAL(move_updates%bias_trans%attempts,dp)/0.5E0_dp
               IF (test_ratio .GT. 2.0E0_dp) test_ratio=2.0E0_dp
               IF (test_ratio .LT. 0.5E0_dp) test_ratio=0.5E0_dp
               rmtrans(molecule_type)=rmtrans(molecule_type)*test_ratio
            ENDIF

 ! make an upper bound...10 a.u.
            IF (rmtrans(molecule_type) .GT. 10.0E0_dp) &
                 rmtrans(molecule_type )= 10.0E0_dp

 ! clear the counters
            move_updates%bias_trans%attempts=0
            move_updates%bias_trans%successes=0

! write the new displacement to a file
            IF(ionode) WRITE(rm,*) nnstep,' rmtrans = ',&
                 rmtrans(molecule_type)*angstrom,' angstroms'
         ENDIF

! update the maximum displacement for molecule rotation
         IF( move_updates%bias_rot%attempts .GT. 0 ) THEN

! first account for the extreme cases
            IF ( move_updates%bias_rot%successes == 0 ) THEN
               rmrot=rmrot/2.0E0_dp

               IF (rmrot(molecule_type) .GT. pi) rmrot(molecule_type)=pi

            ELSEIF ( move_updates%bias_rot%successes == &
                 move_updates%bias_rot%attempts ) THEN
               rmrot(molecule_type)=rmrot(molecule_type)*2.0E0_dp

! more than pi rotation is meaningless
               IF (rmrot(molecule_type) .GT. pi) rmrot(molecule_type)=pi

            ELSE
! now for the middle case
               test_ratio=REAL(move_updates%bias_rot%successes,dp)&
                    /REAL(move_updates%bias_rot%attempts,dp)/0.5E0_dp
               IF (test_ratio .GT. 2.0E0_dp) test_ratio=2.0E0_dp
               IF (test_ratio .LT. 0.5E0_dp) test_ratio=0.5E0_dp
               rmrot(molecule_type)=rmrot(molecule_type)*test_ratio

! more than pi rotation is meaningless
               IF (rmrot(molecule_type) .GT. pi) rmrot(molecule_type)=pi

            ENDIF

! clear the counters
            move_updates%bias_rot%attempts=0
            move_updates%bias_rot%successes=0

! write the new displacement to a file
            IF(ionode) WRITE(rm,*) nnstep,' rmrot = ',&
                 rmrot(molecule_type)/pi*180.0E0_dp,' degrees'
         ENDIF

      CASE ("volume")

! update the maximum displacement for volume displacement
         IF ( move_updates%volume%attempts .NE. 0) THEN

! first account for the extreme cases
            IF ( move_updates%volume%successes == 0 ) THEN
               rmvolume=rmvolume/2.0E0_dp

            ELSEIF ( move_updates%volume%successes == &
                  move_updates%volume%attempts ) THEN
               rmvolume=rmvolume*2.0E0_dp
            ELSE
! now for the middle case
               test_ratio=REAL(move_updates%volume%successes,dp)/&
                  REAL(move_updates%volume%attempts,dp)/0.5E0_dp
               IF (test_ratio .GT. 2.0E0_dp) test_ratio=2.0E0_dp
               IF (test_ratio .LT. 0.5E0_dp) test_ratio=0.5E0_dp
               rmvolume=rmvolume*test_ratio

            ENDIF

! clear the counters
            move_updates%volume%attempts=0
            move_updates%volume%successes=0

! write the new displacement to a file
            IF(ionode) WRITE(rm,*) nnstep,' rmvolume = ',&
                    rmvolume*angstrom**3,' angstroms^3'

         ENDIF

      END SELECT

! set some of the MC parameters
      CALL set_mc_par(mc_par,rmbond=rmbond,rmangle=rmangle,rmrot=rmrot,&
         rmtrans=rmtrans,rmvolume=rmvolume,rmdihedral=rmdihedral)

! end the timing
      CALL timestop(handle)

  END SUBROUTINE mc_move_update

END MODULE mc_move_control

