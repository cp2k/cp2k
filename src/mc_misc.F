!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/mc_misc [1.0] *
!!
!!   NAME
!!     mc_misc
!!
!!   FUNCTION
!!     contains miscellaneous subroutines used in the Monte Carlo runs
!! 
!!   AUTHOR
!!     MJM
!!
!******************************************************************************

MODULE mc_misc
  USE cp_error_handling,               ONLY: cp_error_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE kinds,                           ONLY: dbl,&
                                             dp
  USE mathconstants,                   ONLY: pi
  USE mc_control,                      ONLY: mc_simulation_parameters_type
  USE mc_move_control,                 ONLY: mc_moves_type
  USE particle_list_types,             ONLY: particle_list_type
  USE physcon,                         ONLY: angstrom,&
                                             n_avogadro
  USE qs_parser,                       ONLY: close_file,&
                                             open_file
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE :: nonbonded_energy
  PUBLIC :: get_water_data,final_mc_write,mc_classical_energy,&
            mc_averages_type,mc_averages_p_type,init_mc_averages,&
            get_center_of_mass,mc_make_dat_file,mc_coordinate_fold

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_misc"

  TYPE mc_averages_type
      REAL(KIND = dp) :: ave_energy
      REAL(KIND = dp) :: ave_energy_squared
      REAL(KIND = dp) :: ave_cell_length
      REAL(KIND = dp) :: bond1
      REAL(KIND = dp) :: bond2
      REAL(KIND = dp) :: angle
      REAL(KIND = dp) :: molecules
  END TYPE mc_averages_type
  TYPE mc_averages_p_type
      TYPE (mc_averages_type),POINTER :: averages
  END TYPE mc_averages_p_type

!******************************************************************************

CONTAINS

!******************************************************************************
!!****** mc_misc.F/init_mc_averages*
!!
!!   NAME
!!     init_mc_averages
!!
!!   SYNOPSIS
!!     Subroutine init_mc_averages ( averages )
!!       Type(mc_averages_type), Pointer :: averages
!!     End Subroutine init_mc_averages
!!
!!   FUNCTION
!!     initializes the structure that holds running averages of MC variables
!!
!!   ARGUMENTS
!!     - averages: the mc_averages strucutre you want to initialize
!!
!!   AUTHOR
!!     MJM
!!
!!****************************************************************************

SUBROUTINE init_mc_averages ( averages  )


    TYPE(mc_averages_type), POINTER          :: averages

    CHARACTER(LEN=*), PARAMETER      :: routineN = "init_mc_averages"

    INTEGER                                  :: handle, stat

!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! allocate all the structures...not sure why, but it won't work otherwise
      ALLOCATE (averages,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "averages",0)

      averages%ave_energy=0.0d0
      averages%ave_energy_squared=0.0d0
      averages%ave_cell_length=0.0d0
      averages%bond1=0.0d0
      averages%bond2=0.0d0
      averages%angle=0.0d0
      averages%molecules=0.0D0

! end the timing
       CALL timestop(0.0d0,handle)


END SUBROUTINE init_mc_averages

!******************************************************************************
!!****** mc_misc.F/get_water_data *
!!
!!   NAME
!!     get_water_data
!!
!!   SYNOPSIS
!!     Subroutine get_water_data (force_env,averages,nnstep,nchain,data)
!!       Type(force_env_type), Pointer            :: force_env
!!       Type(mc_averages_type), Pointer          :: averages
!!       Integer, Intent(IN)                      :: nnstep
!!       Integer, Intent(IN)                      :: nchain
!!       Integer, Intent(IN)                      :: data
!!     End Subroutine get_water_data
!!
!!   FUNCTION
!!     gets the running average of both bond lengths and the bond
!!     angle for every water molecule in the system
!!
!!   ARGUMENTS
!!     - force_env: the force_env you want to get information on
!!     - averages: the mc_averages structure to store the information
!!     - nnstep: the number of times this these averages have been
!!               done (including this one)
!!     - nchain: the number of molecules in this force_env
!!     - data: the unit to write the data to
!!
!!   AUTHOR
!!     MJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!*** ************************************************************************

SUBROUTINE get_water_data ( force_env,averages,nnstep,nchain,data  )


    TYPE(force_env_type), POINTER            :: force_env
    TYPE(mc_averages_type), POINTER          :: averages
    INTEGER, INTENT(IN)                      :: nnstep, nchain, data

    CHARACTER(LEN=*), PARAMETER        :: routineN = "get_water_data"

    INTEGER                                  :: handle, iimol, imol
    REAL(KIND=dp)                            :: angle_loc, bond1_loc, &
                                                bond2_loc, tot_squared
    REAL(KIND=dp), DIMENSION(1:3)            :: vec_1, vec_2
    TYPE(cp_error_type)                      :: error
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: oldsys
    TYPE(particle_list_type), POINTER        :: particles

!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! initialize some stuff
      NULLIFY(oldsys,particles)
      bond1_loc=0.0d0
      bond2_loc=0.0d0
      angle_loc=0.0d0

! get the particle coordinates
      CALL force_env_get(force_env,subsys=oldsys,error=error)
      CALL cp_subsys_get(oldsys(1)%subsys, &
            particles=particles, error=error)

! find the running average bond lengths
       DO iimol=1,nchain
          imol=3*iimol-2

          tot_squared=(particles%els(imol)%r(1)-particles%els(imol+1)%&
                 r(1))**2+(particles%els(imol)%r(2)-particles%els(imol&
                 +1)%r(2))**2+(particles%els(imol)%r(3)-particles%els(&
                 imol+1)%r(3))**2
          bond1_loc=bond1_loc*REAL(iimol-1,dbl)/REAL(iimol,dbl) &
                  +DSQRT(tot_squared)/REAL(iimol,dbl)

          tot_squared=(particles%els(imol)%r(1)-particles%els(imol+2)%&
                 r(1))**2+(particles%els(imol)%r(2)-particles%els(imol&
                 +2)%r(2))**2+(particles%els(imol)%r(3)-particles%els(&
                 imol+2)%r(3))**2
          bond2_loc=bond2_loc*REAL(iimol-1,dbl)/REAL(iimol,dbl) &
                  +DSQRT(tot_squared)/REAL(iimol,dbl)

! now find the running average bond angle
          vec_1(1:3)=particles%els(imol+1)%r(1:3)-&
              particles%els(imol)%r(1:3)
          vec_2(1:3)=particles%els(imol+2)%r(1:3)-&
              particles%els(imol)%r(1:3)

          angle_loc=angle_loc*REAL(iimol-1,dbl)/REAL(iimol,dbl) + &
                   1.0d0/REAL(iimol,dbl)*&
              DACOS(DOT_PRODUCT(vec_1,vec_2)/(bond1_loc*bond2_loc))&
              *180.0d0/pi

       ENDDO
       averages%bond1=averages%bond1*REAL(nnstep-1,dbl)/&
                  REAL(nnstep,dbl)+bond1_loc/REAL(nnstep,dbl)
       averages%bond2=averages%bond2*REAL(nnstep-1,dbl)/&
                  REAL(nnstep,dbl)+bond2_loc/REAL(nnstep,dbl)
       averages%angle=averages%angle*REAL(nnstep-1,dbl)/&
                  REAL(nnstep,dbl)+angle_loc/REAL(nnstep,dbl)

       WRITE(data,*) nnstep,bond1_loc*angstrom,bond2_loc*angstrom,&
                     angle_loc
 
! end the timing
       CALL timestop(0.0d0,handle)


END SUBROUTINE get_water_data

!******************************************************************************

!******************************************************************************
!!****s* mc_misc.F/final_mc_write
!!
!!   NAME
!!     final_mc_write
!!
!!   SYNOPSIS
!!     Subroutine final_mc_write (mc_par,moves,iw,energy_check,initial_energy,
!!                           final_energy,averages,move_updates)
!!           Type(mc_simulation_parameters_type),Pointer
!!                                                 :: mc_par
!!           Type(mc_moves_type),Pointer           :: moves
!!           Type(mc_moves_type),Pointer           :: move_updates
!!           Type(mc_averages_type),Pointer        :: averages
!!           Real(Kind=dp), Intent(IN)             :: energy_check
!!           Real(Kind=dp), Intent(IN)             :: initial_energy
!!           Real(Kind=dp), Intent(IN)             :: final_energy
!!           Integer, Intent(IN)                   :: iw
!!     End Subroutine final_mc_write
!!
!!   FUNCTION
!!     writes a bunch of simulation data to the specified unit
!!
!!   ARGUMENTS
!!     - mc_par: the mc parameters for the simulation
!!     - moves: the structure that holds data on how many moves are
!!              accepted/rejected
!!     - move_updates: the structure that holds data on how many classical 
!!              moves are accepted/rejected, used for updating maximum
!!              displacements
!!     - averages: the structure that holds computed average properites for
!!              the simulation
!!     - initial_energy: the initial energy of the system
!!     - final_energy: the final energy of the system
!!     - energy_check: the sum of the energy changes of each move
!!     - iw: the unit to write to
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

SUBROUTINE final_mc_write (mc_par,moves,iw,energy_check,initial_energy,&
                           final_energy,averages,move_updates)


    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(mc_moves_type), POINTER             :: moves,move_updates
    INTEGER, INTENT(IN)                      :: iw
    REAL(KIND=dp), INTENT(IN)         :: energy_check, initial_energy, &
                                                final_energy
    TYPE(mc_averages_type), POINTER          :: averages

    CHARACTER(LEN=*), PARAMETER         :: routineN = "final_mc_write"

    INTEGER                                  :: handle

!------------------------------------------------------------------------------
! begin the timing of the subroutine

    CALL timeset(routineN,'I','',handle)

! write out some data averaged over the whole simulation
    WRITE(iw,*)
    WRITE(iw,*) '****************************************************',&
                   '***************************'
    WRITE(iw,'(A,T66,F15.8)') " Average Energy [Hartrees]:", &
                averages%ave_energy
    WRITE(iw,'(A,T66,F15.8)') " Average O-H distance 1 [angstroms]:",&
                   averages%bond1*angstrom
    WRITE(iw,'(A,T66,F15.8)') " Average O-H distance 2 [angstroms]:",&
                   averages%bond2*angstrom
    WRITE(iw,'(A,T66,F15.8)') " Average H-O-H angle [degrees]:",&
                   averages%angle
    IF (mc_par%pmswap .GT. 0.0D0) &
      WRITE(iw,'(A,T66,F15.8)') " Average number of molecules:",&
                   averages%molecules
    IF (averages%ave_cell_length .NE. 0.0d0) THEN
       WRITE(iw,'(A,A,T68,F13.10)') " Average Density ",&
               "[g/cm^3]:",1.0d0/(averages%ave_cell_length)**3&
                      *mc_par%nchain/&
                      (angstrom)**3*(1.0d8)**3/n_avogadro*18.0184d0
       WRITE(iw,'(A,A,T67,F14.10)') " Average Cell Length ",&
               "[angstroms]:",averages%ave_cell_length*angstrom
    ENDIF
    WRITE(iw,*)

! write out acceptance rates for the moves
    WRITE(iw,*) "Move Data"
    WRITE(iw,*) '----------------------------------------------------',&
                   '---------------------------'
    IF (moves%volume%attempts .GT. 0) THEN
       WRITE(iw,'(A,T43,A)') " Volume Moves",&
          "Attempted       Accepted       Percent"
       WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%volume%attempts,&
             moves%volume%successes,&
             real(moves%volume%successes,dbl)/ &
             real(moves%volume%attempts,dbl)*100.0d0
       WRITE(iw,'(A,T71,F10.5)') &
           ' Maximum Volume displacement [angstroms^3] = ',&
                     mc_par%rmvolume*(angstrom)**3
       WRITE(iw,*) '-----------------------------------------------',&
      '--------------------------------'
    ENDIF

    IF (moves%angle%attempts .GT. 0 .OR. &
                            moves%bond%attempts .GT. 0) THEN
       WRITE(iw,'(A,T43,A)') " Conformational Moves",&
          "Attempted       Accepted       Percent"
       WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
          moves%bond%attempts+moves%angle%attempts,&
          moves%bond%successes+moves%angle%successes,&
          real(moves%bond%successes+moves%angle%successes,dbl)/ &
          real(moves%bond%attempts+moves%angle%attempts,dbl)*100.0d0
       IF (moves%bond%attempts .GT. 0 ) THEN
          WRITE(iw,'(A,T43,A)') " Bond Changes",&
            "Attempted       Accepted       Percent"
          WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
            moves%bond%attempts,&
            moves%bond%successes,&
            real(moves%bond%successes,dbl)/ &
            real(moves%bond%attempts,dbl)*100.0d0
          IF( .NOT. mc_par%lclassical) WRITE(iw,'(A,T71,F10.5)') &
            ' Maximum bond displacement [angstroms]= ',&
            mc_par%rmbond*angstrom
       ENDIF
       IF (moves%angle%attempts .GT. 0 ) THEN
          WRITE(iw,'(A,T43,A)') " Angle Changes",&
            "Attempted       Accepted       Percent"
          WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
            moves%angle%attempts,&
            moves%angle%successes,&
            real(moves%angle%successes,dbl)/ &
            real(moves%angle%attempts,dbl)*100.0d0
          IF( .NOT. mc_par%lclassical) WRITE(iw,'(A,T71,F10.5)') &
            ' Maximum angle displacement [degrees]= ',&
            mc_par%rmangle/pi*180.0d0
       ENDIF
       WRITE(iw,*) '-----------------------------------------------',&
      '--------------------------------'
    ENDIF

    IF (moves%trans%attempts .GT. 0) THEN
       WRITE(iw,'(A,T43,A)') " Translation Moves",&
          "Attempted       Accepted       Percent"
       write(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%trans%attempts,&
             moves%trans%successes,&
             real(moves%trans%successes,dbl)/ &
             real(moves%trans%attempts,dbl)*100.0d0
       IF( .NOT. mc_par%lclassical) WRITE(iw,'(A,T71,F10.5)') &
         ' Maximum molecular translational displacement [angstroms]= ',&
                  mc_par%rmtrans*angstrom
!       WRITE(iw,'(A,T71,F10.5)') &
!          ' Total translational displacement [angstroms]=',&
!                  moves%trans_dis*angstrom
       WRITE(iw,*) '-----------------------------------------------',&
      '--------------------------------'
    ENDIF

    IF (moves%rot%attempts .GT. 0) THEN
       WRITE(iw,'(A,T43,A)') " Rotation Moves",&
          "Attempted       Accepted       Percent"
       WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%rot%attempts,&
             moves%rot%successes,&
             real(moves%rot%successes,dbl)/ &
             real(moves%rot%attempts,dbl)*100.0d0
       IF( .NOT. mc_par%lclassical) WRITE(iw,'(A,T71,F10.5)') &
           ' Maximum rotational displacement [degrees] = ',&
                           mc_par%rmrot/pi*180.0d0
       WRITE(iw,*) '-----------------------------------------------',&
      '--------------------------------'
    ENDIF

    IF (moves%Quickstep%attempts .GT. 0) THEN
       WRITE(iw,'(A,T43,A)') " Quickstep Moves",&
          "Attempted       Accepted       Percent"
       WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%Quickstep%attempts,&
             moves%Quickstep%successes,&
             real(moves%Quickstep%successes,dbl)/ &
             real(moves%Quickstep%attempts,dbl)*100.0d0
       WRITE(iw,*) '-----------------------------------------------',&
      '--------------------------------'
    ENDIF

    IF (moves%swap%attempts .GT. 0) THEN
       WRITE(iw,'(A,T43,A)') " Swap Moves into this box",&
          "Attempted       Empty          Percent"
       WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%swap%attempts,&
             moves%empty,&
             real(moves%empty,dbl)/ &
             real(moves%swap%attempts,dbl)*100.0d0
       WRITE(iw,'(A,T43,A)') "                  Growths",&
          "Attempted       Sucessful      Percent"
       WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%swap%attempts,&
             moves%grown,&
             real(moves%grown,dbl)/ &
             real(moves%swap%attempts,dbl)*100.0d0
       WRITE(iw,'(A,T43,A)') "                    Total",&
          "Attempted       Accepted       Percent"
       WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%swap%attempts,&
             moves%swap%successes,&
             real(moves%swap%successes,dbl)/ &
             real(moves%swap%attempts,dbl)*100.0d0
       WRITE(iw,*) '-----------------------------------------------',&
      '--------------------------------'
    ENDIF

! now we write out information on the classical moves, if it's
! a classical simulations
    IF (mc_par%lclassical) THEN
       WRITE(iw,*) "Classical Move Data"
       WRITE(iw,*) '------------------------------------------------',&
                   '-------------------------------'
! first, bond changes...record all the moves done since the last update,
! and print out the totals
       IF (moves%bond%attempts .GT. 0 ) THEN
          move_updates%classical_bond%attempts=&
             move_updates%classical_bond%attempts+&
             move_updates%bond%attempts
          move_updates%classical_bond%successes=&
             move_updates%classical_bond%successes+&
             move_updates%bond%successes
          WRITE(iw,'(A,T43,A)') " Bond Changes",&
            "Attempted       Accepted       Percent"
          WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
            move_updates%classical_bond%attempts,&
            move_updates%classical_bond%successes,&
            real(move_updates%classical_bond%successes,dbl)/ &
            real(move_updates%classical_bond%attempts,dbl)*100.0d0
          WRITE(iw,'(A,T71,F10.5)') &
            ' Maximum bond displacement [angstroms]= ',&
            mc_par%rmbond*angstrom
          WRITE(iw,*) '-------------------------------------------',&
                   '------------------------------------'
       ENDIF

! now angle changes
       IF (moves%angle%attempts .GT. 0 ) THEN
          move_updates%classical_angle%attempts=&
             move_updates%classical_angle%attempts+&
             move_updates%angle%attempts
          move_updates%classical_angle%successes=&
             move_updates%classical_angle%successes+&
             move_updates%angle%successes
          WRITE(iw,'(A,T43,A)') " Angle Changes",&
            "Attempted       Accepted       Percent"
          WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
            move_updates%classical_angle%attempts,&
            move_updates%classical_angle%successes,&
            real(move_updates%classical_angle%successes,dbl)/ &
            real(move_updates%classical_angle%attempts,dbl)*100.0d0
          WRITE(iw,'(A,T71,F10.5)') &
            ' Maximum angle displacement [degrees]= ',&
            mc_par%rmangle/pi*180.0d0
          WRITE(iw,*) '-------------------------------------------',&
                   '------------------------------------'
       ENDIF

! now molecule translations
       IF (moves%trans%attempts .GT. 0 ) THEN
          move_updates%classical_trans%attempts=&
             move_updates%classical_trans%attempts+&
             move_updates%trans%attempts
          move_updates%classical_trans%successes=&
             move_updates%classical_trans%successes+&
             move_updates%trans%successes
          WRITE(iw,'(A,T43,A)') " Molecule Translations",&
            "Attempted       Accepted       Percent"
          WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
            move_updates%classical_trans%attempts,&
            move_updates%classical_trans%successes,&
            real(move_updates%classical_trans%successes,dbl)/ &
            real(move_updates%classical_trans%attempts,dbl)*100.0d0
          WRITE(iw,'(A,T71,F10.5)') &
         ' Maximum molecular translational displacement [angstroms]= ',&
                  mc_par%rmtrans*angstrom
          WRITE(iw,*) '-------------------------------------------',&
                   '------------------------------------'
       ENDIF

! now molecule rotations
       IF (moves%rot%attempts .GT. 0 ) THEN
          move_updates%classical_rot%attempts=&
             move_updates%classical_rot%attempts+&
             move_updates%rot%attempts
          move_updates%classical_rot%successes=&
             move_updates%classical_rot%successes+&
             move_updates%rot%successes
          WRITE(iw,'(A,T43,A)') " Molecule Rotations",&
            "Attempted       Accepted       Percent"
          WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
            move_updates%classical_rot%attempts,&
            move_updates%classical_rot%successes,&
            real(move_updates%classical_rot%successes,dbl)/ &
            real(move_updates%classical_rot%attempts,dbl)*100.0d0
          WRITE(iw,'(A,T71,F10.5)') &
           ' Maximum rotational displacement [degrees] = ',&
                           mc_par%rmrot/pi*180.0d0
          WRITE(iw,*) '-------------------------------------------',&
                   '------------------------------------'
       ENDIF

    ENDIF



! see if the energies add up properly
    IF(ABS(initial_energy+energy_check-final_energy) .GT. 0.0000001d0) &
         THEN
         WRITE(iw,*) '!!!!!!! We have an energy problem. !!!!!!!!'
         WRITE(iw,'(A,T64,F16.10)') 'Final Energy = ',final_energy
         WRITE(iw,'(A,T64,F16.10)') 'Inital Energy + energy_check =',&
                 initial_energy+energy_check
    ENDIF
    WRITE(iw,*) '****************************************************',&
                   '***************************'
    WRITE(iw,*)
 
! end the timing
    CALL timestop(0.0d0,handle)

 
END SUBROUTINE final_mc_write

!******************************************************************************
!!****s* mc_misc.F/mc_classical_energy
!!
!!   NAME
!!     mc_classical_energy
!!
!!   SYNOPSIS
!!     Subroutine mc_classical_energy ( force_env, classical_energy_new,
!!                                      mc_par, nstart, nend, rcut)
!!           Type(mc_simulation_parameters_type),Pointer
!!                                                 :: mc_par
!!           Type(force_env_type), Pointer         :: force_env
!!           Real(Kind=dp), Intent(IN)             :: rcut
!!           Real(Kind=dp), Intent(OUT)            :: classical_energy_new
!!           Integer, Intent(IN)                   :: nstart
!!           Integer, Intent(IN)                   :: nend
!!     End Subroutine mc_classical_energy
!!
!!   FUNCTION
!!     calculates the classical energy of a system configuration using
!!     a potential force-matched to ab initio results for water
!!
!!   ARGUMENTS
!!     - mc_par: the mc parameters for the force env
!!     - force_env: the force environment containing the coordinates
!!     - rcut: the intermolecular cutoff you want to use for the potential,
!!             in atomic units...be careful with this, as the potential is
!!             only zero at the cutoff for certain values
!!     - nstart: the number of the molecule start calculating for
!!     - nend: the number of the moelcule to end calculating for
!!     --- if nstart is not equal to 1 and nend is not equal to the total
!!         number of molecules, set them equal
!!     - classical_energy_new: the classical energy returned by the routine
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE mc_classical_energy ( force_env, classical_energy_new,&
                                 mc_par, nstart, nend ,rcut)

    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), INTENT(OUT)               :: classical_energy_new
    REAL(KIND=dp), INTENT(IN)               :: rcut
    INTEGER, INTENT(IN)               :: nstart,nend
    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_classical_energy"

    INTEGER                                  :: istart,iend,jstart,jend
    INTEGER                                  :: handle,imol,ipart,jmol,&
                                                jpart,ntij,i,o,m,n
    REAL(KIND=dp) :: box_length, charge_energy, dist, &
              LJ_energy, unit_energy,&
              stretch_energy,bend_energy,mag_HH,kstretch,OH_bond1,&
              OH_bond2
    REAL(KIND=dp), DIMENSION(1:3)            :: charges, min_IJ, RIJ,&
                                 kbend,HH,OH1,OH2
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: oldsys
    TYPE(cp_error_type)                      :: error
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), &
      DIMENSION(1:3, 1:3, 1:mc_par%nchain)   :: r
    TYPE(particle_list_type), POINTER        :: particles
    LOGICAL                                  :: lcalc_charge,lall

!------------------------------------------------------------------------------
! begin the timing of the subroutine
      CALL timeset(routineN,'I','',handle)


! initialize some stuff
      NULLIFY(oldsys,particles)
      LJ_energy=0.0d0
      charge_energy=0.0d0
      bend_energy=0.0d0
      stretch_energy=0.0d0
      kstretch= 0.4015204d0 ! Hartrees/a.u.**2
      kbend(1)= 0.02384043d0 ! Hartrees/a.u.**2
      kbend(2)= -0.0506033d0 ! Hartrees/a.u.**3
      kbend(3)= 0.0287044d0  ! Hartrees/a.u.**4
      charges(1:3) = (/ -0.7699784d0,0.3849892d0,0.3849892d0 /)
!      charges(:)=charges(:)*1.60d-19 ! Coulumbs

! get the particle coordinates and the cell length
      CALL force_env_get(force_env,cell=cell,subsys=oldsys,error=error)
      CALL cp_subsys_get(oldsys(1)%subsys, &
            particles=particles, error=error)

      box_length=cell%hmat(1,1)

!      IF (box_length .GT. 15.6d0/angstrom) THEN
!         rcut=7.8d0/angstrom    ! always the same for this potential
!      ELSE
!         rcut=box_length/2.0d0
!      ENDIF

!      rcut=7.8d0/angstrom ! we're now using multiple images

! put the coordinates into an easier matrix to manipulate
      DO imol=1,mc_par%nchain
         r(1:3,1,imol)=particles%els(imol*3-2)%r(1:3)
         r(1:3,2,imol)=particles%els(imol*3-1)%r(1:3)
         r(1:3,3,imol)=particles%els(imol*3)%r(1:3)
!         if(mc_par%ionode) WRITE(6,*) 'One ',r(1:3,1,imol)*angstrom
!         if(mc_par%ionode) WRITE(6,*) 'Two ',r(1:3,2,imol)*angstrom
!         if(mc_par%ionode) WRITE(6,*) 'Three ',r(1:3,3,imol)*angstrom
      ENDDO

! we need to calculate the energy of this configuration
! first, let's calculate intramolecular (bending and stretching)
      DO imol=nstart,nend
         
! stretching
         OH_bond1=DSQRT((r(1,1,imol)-r(1,2,imol))**2+&
                 (r(2,1,imol)-r(2,2,imol))**2+&
                   (r(3,1,imol)-r(3,2,imol))**2)
         OH_bond2=DSQRT((r(1,1,imol)-r(1,3,imol))**2+&
                   (r(2,1,imol)-r(2,3,imol))**2+&
                   (r(3,1,imol)-r(3,3,imol))**2)

         stretch_energy=stretch_energy+0.5d0*kstretch*((OH_bond1-&
            1.81200d0)**2+(OH_bond2-1.81200d0)**2)

! bond bending
         DO i=1,3
            HH(i)=r(i,2,imol)-r(i,3,imol)
         ENDDO
         mag_HH=DSQRT(DOT_PRODUCT(HH,HH))
!         bend_energy=bend_energy+0.5d0*kbend(1)*(mag_HH-2.88799d0)**2+&
!             0.33333333333333333d0*kbend(2)*(mag_HH-2.88799d0)**3+&
!             0.25d0*kbend(3)*(mag_HH-2.88799d0)**4
         bend_energy=bend_energy+0.5d0*kbend(1)*(mag_HH-3.0699d0)**2+&
             0.33333333333333333d0*kbend(2)*(mag_HH-3.0699d0)**3+&
             0.25d0*kbend(3)*(mag_HH-3.0699d0)**4
!         bend_energy=bend_energy+0.5d0*kbend(1)*(mag_HH-3.36d0)**2+&
!             0.25d0*kbend(3)*(mag_HH-3.36d0)**4
!         bend_energy=bend_energy*2.0d0

      ENDDO
!      bend_energy=bend_energy*4.0d0

! now let's find the LJ energy between all the oxygens and
! the charge interactions
      lall=.FALSE.
      jstart=1
      IF( nstart == 1 .AND. nend == mc_par%nchain) lall=.TRUE.
      DO imol=nstart,nend
         IF (lall) jstart=imol+1
         DO jmol=jstart,mc_par%nchain
            IF(imol == jmol) CYCLE

! now let's triple loop over all the boxes touching the simulation box
!            DO m=-1,1
!                DO n=-1,1
!                  DO o=-1,1

! find the minimum image distance between oxygens
            RIJ(1)=r(1,1,imol)-r(1,1,jmol)-box_length*ANINT(&
               (r(1,1,imol)-r(1,1,jmol))/box_length)
            RIJ(2)=r(2,1,imol)-r(2,1,jmol)-box_length*ANINT(&
               (r(2,1,imol)-r(2,1,jmol))/box_length)
            RIJ(3)=r(3,1,imol)-r(3,1,jmol)-box_length*ANINT(&
               (r(3,1,imol)-r(3,1,jmol))/box_length)
!            RIJ(1)=r(1,1,imol)-r(1,1,jmol)-box_length*m
!            RIJ(2)=r(2,1,imol)-r(2,1,jmol)-box_length*n
!            RIJ(3)=r(3,1,imol)-r(3,1,jmol)-box_length*o
            
            dist=DSQRT(DOT_PRODUCT(RIJ,RIJ))

            IF (dist < rcut ) THEN
!               ipart=1
!               jpart=1
               lcalc_charge = .true.
               CALL nonbonded_energy(dist,unit_energy,1)
               LJ_energy=LJ_energy+unit_energy               
               charge_energy=charge_energy+charges(1)&
                   *charges(1)/(dist) 
            ELSE
               lcalc_charge=.false.
            ENDIF

            DO ipart=1,3
               DO jpart=1,3        
                  IF (ipart == 1 .AND. jpart == 1) CYCLE

! find the minimum image distance
                  RIJ(1)=r(1,ipart,imol)-r(1,jpart,jmol)-box_length&
                    *ANINT((r(1,ipart,imol)-r(1,jpart,jmol))/box_length)
                  RIJ(2)=r(2,ipart,imol)-r(2,jpart,jmol)-box_length&
                    *ANINT((r(2,ipart,imol)-r(2,jpart,jmol))/box_length)
                  RIJ(3)=r(3,ipart,imol)-r(3,jpart,jmol)-box_length&
                    *ANINT((r(3,ipart,imol)-r(3,jpart,jmol))/box_length)
!                  RIJ(1)=r(1,ipart,imol)-r(1,jpart,jmol)-box_length*m
!                  RIJ(2)=r(2,ipart,imol)-r(2,jpart,jmol)-box_length*n
!                  RIJ(3)=r(3,ipart,imol)-r(3,jpart,jmol)-box_length*o                  
                  dist=DSQRT(DOT_PRODUCT(RIJ,RIJ))

! LJ energy
                  IF ( dist < rcut ) THEN
                     ntij=(ipart-1)*10 + jpart
                     CALL nonbonded_energy(dist,unit_energy,ntij)
                     LJ_energy=LJ_energy+unit_energy
                  ENDIF

! coulombic
                  IF ( lcalc_charge ) THEN
                     charge_energy=charge_energy&
                        +charges(ipart)*charges(jpart)/(dist)
                  ENDIF

               ENDDO
            ENDDO
 !        ENDDO
 !     ENDDO
 !     ENDDO
         ENDDO
      ENDDO

! sum up the total energy
      classical_energy_new=bend_energy+&
                stretch_energy+LJ_energy+charge_energy
                
!      WRITE(mc_par%diff,*) 'bending',bend_energy,classical_energy_new 
!      classical_energy_new=LJ_energy+charge_energy
!      IF (nstart == 12) THEN
!      if(mc_par%ionode) WRITE(6,*) 'enegies =',classical_energy_new/3.16682968D-6,LJ_energy/3.16682968D-6,&
!                 charge_energy/3.16682968D-6,bend_energy/3.16682968D-6,&
!                 stretch_energy/3.16682968D-6
!      STOP
!      ENDIF
! scale the total energy to get classical vs. Quickstep energy diffs closer
!      classical_energy_new=classical_energy_new*1.5d0

! end the timing
      CALL timestop(0.0d0,handle)


END SUBROUTINE mc_classical_energy

!******************************************************************************
!!****s* mc_misc.F/nonbonded_energy
!!
!!   NAME
!!     nonbonded_energy
!!
!!   SYNOPSIS
!!     Subroutine nonbonded_energy(dist,ljfake,ntij)
!!           Real(Kind=dp), Intent(IN)             :: dist
!!           Real(Kind=dp), Intent(OUT)            :: ljfake
!!           Integer, Intent(IN)                   :: ntij
!!     End Subroutine nonbonded_energy
!!
!!   FUNCTION
!!     calculates the nonbonded energy of two particles, based on the 
!!     seperation and atom types
!!
!!   ARGUMENTS
!!     - dist: the distance between the atoms (in a.u.)
!!     - ljfake: the energy generated by the subroutine (fake Lennard-Jones 
!!               potential), in Hartree
!!     - ntij: a number that describes the interaction
!!             1: oxygen-oxygen
!!             2,3,11,21: oxygen-hydrogen
!!             12,13,22,23L hydrogen-hydrogen
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE nonbonded_energy(dist,ljfake,ntij)

      IMPLICIT NONE
      INTEGER,INTENT(IN) :: ntij
      REAL(KIND=dp),INTENT(IN) :: dist
      REAL(KIND=dp),DIMENSION(2:15,1:3) :: A
      REAL(KIND=dp),DIMENSION(1:3) :: r_core,pshift,b,m
      REAL(KIND=dp),INTENT(OUT) :: ljfake
      INTEGER :: i,j,index
      
      index=0
      IF (ntij == 1 ) index=1
      IF (ntij == 2 ) index=2
      IF (ntij == 3 ) index=2
      IF (ntij == 11 ) index=2
      IF (ntij == 12 ) index=3
      IF (ntij == 13 ) index=3
      IF (ntij == 21 ) index=2
      IF (ntij == 22 ) index=3
      IF (ntij == 23 ) index=3
      IF (index == 0) WRITE(6,*) 'index =',index

! --- these pshift make the potential zero at 7.8 angstroms
!      pshift(1)=0.127827d0 ! Hartree
!      pshift(2)=2.22408d0
!      pshift(3)=0.256530d0
! --- these make it zero at 6.211 angstroms
      pshift(1)=0.128426d0 ! Hartree
      pshift(2)=2.223850d0
      pshift(3)=0.256652d0
      m(1) = 0.0d0   ! Hartree/a.u.
      m(2) = -0.00437054d0
      m(3) = 0.00377391d0
! --- these b make the potential zero at 7.8 angstroms
!      b(1) = 1.0d23 ! Hartree
!      b(2) = 0.0313652d0
!      b(3) = -0.0206204d0
! --- these b make the potential zero at 6.211 angstroms
      b(1) = 1.0d23 ! Hartree
      b(2) = 0.0311350d0
      b(3) = -0.0204983d0
      r_core(1)=0.9d0 ! a.u.
      r_core(2)=2.95d0
      r_core(3)=3.165d0
      A(2,1) = 4.786380682394d0   ! Hartree*a.u.^2
      A(2,2) = -195.7716013277d0
      A(2,3) = -26.29456010782d0
      A(3,1) = -1543.407053545d0
      A(3,2) = 15343.78613395d0
      A(3,3) = 2373.352548248d0
      A(4,1) = 88783.31188529d0
      A(4,2) = -530864.4586516d0
      A(4,3) = -93880.43551360d0
      A(5,1) = -2361200.155376d0
      A(5,2) = 10707934.39058d0
      A(5,3) = 2154624.884809d0
      A(6,1) = 35940504.84679d0
      A(6,2) = -140099704.7890d0
      A(6,3) = -31965151.34955d0
      A(7,1) = -339762743.6358d0
      A(7,2) = 1250943273.785d0
      A(7,3) = 322781785.3278d0
      A(8,1) = 2043874926.466d0
      A(8,2) = -7795458330.676d0
      A(8,3) = -2271097368.668d0
      A(9,1) = -7654856796.383d0
      A(9,2) = 33955897217.31d0
      A(9,3) = 11169163192.90d0
      A(10,1) = 16195251405.65d0
      A(10,2) = -101135640744.0d0
      A(10,3) = -37684457778.47d0
      A(11,1) = -13140392992.18d0
      A(11,2) = 193107995718.7d0
      A(11,3) = 82562104256.03d0
      A(12,1) = -9285572894.245d0
      A(12,2) = -193440560940.0d0
      A(12,3) = -100510435213.4d0
      A(13,1) = 8756947519.029d0
      A(13,2) = -4224406093.918d0
      A(13,3) = 24570342714.65d0
      A(14,1) = 15793297761.67d0
      A(14,2) = 217192386506.5d0
      A(14,3) = 88766181532.94d0
      A(15,1) = 12917180227.21d0
      A(15,2) = -157581228915.5d0
      A(15,3) = -79705131323.98d0
      
      ljfake=0.0d0
      IF (dist .GT. r_core(index)) THEN
! use a series expansion of 1/r
         DO j=2,15
            ljfake=ljfake+A(j,index)/(dist**(j-1)*REAL(j-1,dbl))
        ENDDO
         ljfake=ljfake+pshift(index)
      ELSE
! use a linear potential
         ljfake=m(index)*dist+b(index)
      ENDIF

END SUBROUTINE nonbonded_energy

!******************************************************************************
!!****s* mc_misc.F/get_center_of_mass
!!
!!   NAME
!!     get_center_of_mass
!!
!!   SYNOPSIS
!!     Subroutine get_center_of_mass (coordinates,natom,center_of_mass,mass)
!!           Real(Kind=dp), Dimension(1:3),Intent(OUT)
!!                                                 :: center_of_mass
!!           Real(Kind=dp), Dimension(:,:),Intent(IN)
!!                                                 :: coordinates
!!           Real(Kind=dp), Dimension(:),Intent(IN)
!!                                                 :: mass
!!           Integer, Intent(IN)                   :: natom
!!     End Subroutine get_center_of_mass
!!
!!   FUNCTION
!!     calculates the center of mass of a given molecule
!!
!!   ARGUMENTS
!!     - natom: the number of atoms in the molecule
!!     - coordinates: the coordiantes of the atoms in the molecule
!!     - mass: the mass of the atoms in the molecule
!!     - center_of_mass: the coordinates of the center of mass
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE get_center_of_mass ( coordinates,natom,center_of_mass,&
            mass)

    REAL(KIND=dp), DIMENSION(1:3),INTENT(OUT)  :: center_of_mass
    REAL(KIND=dp), DIMENSION(:,:), INTENT(IN)   :: coordinates
    REAL(KIND=dp), DIMENSION(:),INTENT(IN)     :: mass
    INTEGER, INTENT(IN)               :: natom

    CHARACTER(LEN=*), PARAMETER :: routineN = "get_center_of_mass"

! local variables
    INTEGER                        :: handle,i,iatom
    REAL(KIND=dp)                 :: total_mass
!------------------------------------------------------------------------------
! begin the timing of the subroutine
      CALL timeset(routineN,'I','',handle)

      total_mass=SUM(mass(:))
      center_of_mass(:)=0.0D0

      DO iatom=1,natom
         DO i=1,3
            center_of_mass(i)=center_of_mass(i)+&
               mass(iatom)*coordinates(i,iatom)
         ENDDO
      ENDDO

      center_of_mass(1:3)=center_of_mass(1:3)/total_mass

! end the timing
      CALL timestop(0.0d0,handle)


END SUBROUTINE get_center_of_mass

!******************************************************************************
!!****s* mc_misc.F/get_center_of_mass
!!
!!   NAME
!!     get_center_of_mass
!!
!!   SYNOPSIS
!!     Subroutine get_center_of_mass (coordinates,natom,center_of_mass,mass)
!!           Real(Kind=dp), Dimension(1:3),Intent(OUT)
!!                                                 :: center_of_mass
!!           Real(Kind=dp), Dimension(:,:),Intent(IN)
!!                                                 :: coordinates
!!           Real(Kind=dp), Dimension(:),Intent(IN)
!!                                                 :: mass
!!           Integer, Intent(IN)                   :: natom
!!     End Subroutine get_center_of_mass
!!
!!   FUNCTION
!!     calculates the center of mass of a given molecule
!!
!!   ARGUMENTS
!!     - natom: the total number of atoms
!!     - coordinates: the coordiantes of the atoms in the molecule
!!     - mass: the mass of the atoms in the molecule
!!     - center_of_mass: the coordinates of the center of mass
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE mc_make_dat_file(coordinates,natoms,box_length,filename)

    REAL(KIND=dp), DIMENSION(:,:), INTENT(INOUT)   :: coordinates
    REAL(KIND=dp), INTENT(IN)   :: box_length
    INTEGER, INTENT(IN)                         :: natoms
    CHARACTER(LEN=*),INTENT(IN)                 :: filename

    CHARACTER(LEN=*), PARAMETER :: routineN = "get_center_of_mass"

! local variables
    INTEGER                        :: handle,iatom,unit
    REAL(KIND=dp)                 :: total_mass
!------------------------------------------------------------------------------
! begin the timing of the subroutine
      CALL timeset(routineN,'I','',handle)

! open the file
      CALL open_file(file_name=filename,unit_number=unit,&
         file_action='WRITE',file_status='REPLACE')

! write all the stuff down...this is very specific!
      WRITE(unit,*) "&CP2K"
      WRITE(unit,*) 'PROGRAM      Quickstep'
      WRITE(unit,*) 'PROJECT      H2O_MC'
      WRITE(unit,*) 'IOLEVEL      0'
      WRITE(unit,*) 'FFTLIB       FFTSG'
      WRITE(unit,*) 'RUN_TYPE     MC'
      WRITE(unit,*) "&END"
      WRITE(unit,*) 
      WRITE(unit,*) "&PRINT silent"
      WRITE(unit,*) "&END"
      WRITE(unit,*) 
      WRITE(unit,*) "&MC"
      WRITE(unit,*) "&END"
      WRITE(unit,*) 
      WRITE(unit,*) "&DFT"
      WRITE(unit,*) ' FUNCTIONAL   BLYP'
      WRITE(unit,*) ' XC_SMOOTH_RHO   NN10'
      WRITE(unit,*) ' XC_DERIV      SPLINE2_smooth'
      WRITE(unit,*) "&END"
      WRITE(unit,*) 
      WRITE(unit,*) "&QS"
      WRITE(unit,*) ' CUTOFF       280'
      WRITE(unit,*) "&END"
      WRITE(unit,*) 
      WRITE(unit,*) "&SCF"
      WRITE(unit,*) ' GUESS        ATOMIC'
      WRITE(unit,*) "&END"
      WRITE(unit,*) 
      WRITE(unit,*) "&KIND H"
      WRITE(unit,*) ' BASIS_SET    TZV2P-GTH-BLYP'
      WRITE(unit,*) ' POTENTIAL    GTH-BLYP-q1'
      WRITE(unit,*) "&END"
      WRITE(unit,*)  
      WRITE(unit,*) "&KIND O"
      WRITE(unit,*) ' BASIS_SET    TZV2P-GTH-BLYP'
      WRITE(unit,*) ' POTENTIAL    GTH-BLYP-q6'
      WRITE(unit,*) "&END"
      WRITE(unit,*) 
      WRITE(unit,*) "&CELL"
      WRITE(unit,*) ' UNIT         ANGSTROM'
      WRITE(unit,'(A,3(F7.3,2X))') ' ABC          ',&
         box_length,box_length,box_length
      WRITE(unit,*) "&END"
      WRITE(unit,*) 
      WRITE(unit,*) "&COORD"
      DO iatom=1,natoms,3
         WRITE(unit,'(A,3(F22.15,2X),A)') ' O  ',&
            coordinates(1:3,iatom)*angstrom,'  H2O'
         WRITE(unit,'(A,3(F22.15,2X),A)') ' H  ',&
            coordinates(1:3,iatom+1)*angstrom,'  H2O'
         WRITE(unit,'(A,3(F22.15,2X),A)') ' H  ',&
            coordinates(1:3,iatom+2)*angstrom,'  H2O'
      ENDDO
      WRITE(unit,*) "&END"

! close the file
      CALL close_file(unit_number=unit)

! end the timing
      CALL timestop(0.0d0,handle)


END SUBROUTINE mc_make_dat_file

!******************************************************************************
!!****s* mc_misc.F/get_center_of_mass
!!
!!   NAME
!!     get_center_of_mass
!!
!!   SYNOPSIS
!!     Subroutine get_center_of_mass (coordinates,natom,center_of_mass,mass)
!!           Real(Kind=dp), Dimension(1:3),Intent(OUT)
!!                                                 :: center_of_mass
!!           Real(Kind=dp), Dimension(:,:),Intent(IN)
!!                                                 :: coordinates
!!           Real(Kind=dp), Dimension(:),Intent(IN)
!!                                                 :: mass
!!           Integer, Intent(IN)                   :: natom
!!     End Subroutine get_center_of_mass
!!
!!   FUNCTION
!!     calculates the center of mass of a given molecule
!!
!!   ARGUMENTS
!!     - natom: the total number of atoms
!!     - coordinates: the coordiantes of the atoms in the molecule
!!     - mass: the mass of the atoms in the molecule
!!     - center_of_mass: the coordinates of the center of mass
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE mc_coordinate_fold(coordinates,natoms,box_length)

    REAL(KIND=dp), DIMENSION(:,:), INTENT(INOUT)   :: coordinates
    INTEGER, INTENT(IN)                         :: natoms
    REAL(KIND=dp),INTENT(IN)                 :: box_length

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_coordinate_fold"

! local variables
    INTEGER                        :: handle,iatom,i
    REAL(KIND=dp),DIMENSION(1:3)        :: center_of_mass
    REAL(KIND=dp),DIMENSION(1:3)        :: mass
!------------------------------------------------------------------------------
! begin the timing of the subroutine
      CALL timeset(routineN,'I','',handle)

      mass(1:3)=(/ 15.999d0,1.0074d0,1.0074d0 /)

! loop over all molecules
      DO iatom=1,natoms,3
         CALL get_center_of_mass(coordinates(:,iatom:iatom+2),3,&
            center_of_mass(:),mass(:))
         DO i=1,3
            coordinates(i,iatom)=coordinates(i,iatom)-&
               box_length*FLOOR(center_of_mass(i)/box_length)
            coordinates(i,iatom+1)=coordinates(i,iatom+1)-&
               box_length*FLOOR(center_of_mass(i)/box_length)
            coordinates(i,iatom+2)=coordinates(i,iatom+2)-&
               box_length*FLOOR(center_of_mass(i)/box_length)
         ENDDO

      ENDDO



! end the timing
      CALL timestop(0.0d0,handle)


END SUBROUTINE mc_coordinate_fold

END MODULE mc_misc

!******************************************************************************
