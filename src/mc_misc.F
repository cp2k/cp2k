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
!!     contains miscellaneous subroutines used in the Monte Carlo runs,
!!     mostly I/O stuff
!! 
!!   AUTHOR
!!     MJM
!!
!******************************************************************************

MODULE mc_misc
  USE cp_error_handling,               ONLY: cp_error_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE cp_control_types,                ONLY: dft_control_type
!  USE cp_control_utils,                ONLY: xc_get_rho_smooth_label,&
!                                             xc_get_deriv_method_label
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE kinds,                           ONLY: dp,&
                                             dp
  USE mathconstants,                   ONLY: pi
  USE mc_types,                        ONLY: mc_averages_type,&
                                             mc_moves_type,&
                                             mc_simulation_parameters_type,&
                                             get_mc_par
  USE particle_list_types,             ONLY: particle_list_type
  USE physcon,                         ONLY: angstrom,&
                                             n_avogadro
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_parser,                       ONLY: close_file,&
                                             open_file
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PUBLIC :: final_mc_write,mc_averages_create,mc_averages_release,&
            mc_make_dat_file

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_misc"

!******************************************************************************

CONTAINS

!******************************************************************************
!!****** mc_misc.F/mc_averages_create*
!!
!!   NAME
!!     mc_averages_create
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

SUBROUTINE mc_averages_create ( averages  )


    TYPE(mc_averages_type), POINTER          :: averages

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "mc_averages_create"

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


END SUBROUTINE mc_averages_create

!******************************************************************************
!!****** mc_misc.F/init_mc_averages*
!!
!!   NAME
!!     mc_averages_release
!!
!!   FUNCTION
!!     deallocates the structure that holds running averages of MC variables
!!
!!   ARGUMENTS
!!     - averages: the mc_averages strucutre you want to release
!!
!!   AUTHOR
!!     MJM
!!
!!****************************************************************************

SUBROUTINE mc_averages_release ( averages  )


    TYPE(mc_averages_type), POINTER          :: averages

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "release_mc_averages"

    INTEGER                                  :: handle, stat

!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! deallocate 
      DEALLOCATE (averages,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                     "averages")

      NULLIFY(averages)

! end the timing
       CALL timestop(0.0d0,handle)


END SUBROUTINE mc_averages_release

!******************************************************************************
!!****s* mc_misc.F/final_mc_write
!!
!!   NAME
!!     final_mc_write
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
    TYPE(mc_moves_type), POINTER             :: moves
    INTEGER, INTENT(IN)                      :: iw
    REAL(KIND=dp), INTENT(IN)                :: energy_check, initial_energy, &
                                                final_energy
    TYPE(mc_averages_type), POINTER          :: averages
    TYPE(mc_moves_type), POINTER             :: move_updates

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "final_mc_write"

    LOGICAL                                  :: lclassical
    INTEGER                                  :: handle,nvolmoves
    REAL(KIND=dp)                            :: pmswap,rmvolume,rmbond,&
      rmangle,rmtrans,rmrot
!------------------------------------------------------------------------------
! begin the timing of the subroutine

    CALL timeset(routineN,'I','',handle)

    CALL get_mc_par(mc_par,pmswap=pmswap,rmvolume=rmvolume,&
      lclassical=lclassical,rmbond=rmbond,rmangle=rmangle,&
      rmtrans=rmtrans,rmrot=rmrot,nvolmoves=nvolmoves)

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
    IF (pmswap .GT. 0.0D0) THEN
       WRITE(iw,'(A,T66,F15.8)') " Average number of molecules:",&
                   averages%molecules
    ENDIF
    WRITE(iw,'(A,A,T68,F13.10)') " Average Density ",&
               "[g/cm^3]:",1.0d0/(averages%ave_cell_length)**3&
                      *averages%molecules/&
                      (angstrom)**3*(1.0d8)**3/n_avogadro*18.0184d0
    WRITE(iw,'(A,A,T67,F14.10)') " Average Cell Length ",&
               "[angstroms]:",averages%ave_cell_length*angstrom

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
             REAL(moves%volume%successes,dp)/ &
             REAL(moves%volume%attempts,dp)*100.0d0
       WRITE(iw,'(A,T71,F10.5)') &
           ' Maximum Volume displacement [angstroms^3] = ',&
                     rmvolume*(angstrom)**3
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
          REAL(moves%bond%successes+moves%angle%successes,dp)/ &
          REAL(moves%bond%attempts+moves%angle%attempts,dp)*100.0d0
       IF (moves%bond%attempts .GT. 0 ) THEN
          WRITE(iw,'(A,T43,A)') " Bond Changes",&
            "Attempted       Accepted       Percent"
          WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
            moves%bond%attempts,&
            moves%bond%successes,&
            REAL(moves%bond%successes,dp)/ &
            REAL(moves%bond%attempts,dp)*100.0d0
          IF( .NOT. lclassical) WRITE(iw,'(A,T71,F10.5)') &
            ' Maximum bond displacement [angstroms]= ',&
            rmbond*angstrom
       ENDIF
       IF (moves%angle%attempts .GT. 0 ) THEN
          WRITE(iw,'(A,T43,A)') " Angle Changes",&
            "Attempted       Accepted       Percent"
          WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
            moves%angle%attempts,&
            moves%angle%successes,&
            REAL(moves%angle%successes,dp)/ &
            REAL(moves%angle%attempts,dp)*100.0d0
          IF( .NOT. lclassical) WRITE(iw,'(A,T71,F10.5)') &
            ' Maximum angle displacement [degrees]= ',&
            rmangle/pi*180.0d0
       ENDIF
       WRITE(iw,*) '-----------------------------------------------',&
      '--------------------------------'
    ENDIF

    IF (moves%trans%attempts .GT. 0) THEN
       WRITE(iw,'(A,T43,A)') " Translation Moves",&
          "Attempted       Accepted       Percent"
       WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%trans%attempts,&
             moves%trans%successes,&
             REAL(moves%trans%successes,dp)/ &
             REAL(moves%trans%attempts,dp)*100.0d0
       IF( .NOT. lclassical) WRITE(iw,'(A,T71,F10.5)') &
         ' Maximum molecular translational displacement [angstroms]= ',&
                  rmtrans*angstrom
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
             REAL(moves%rot%successes,dp)/ &
             REAL(moves%rot%attempts,dp)*100.0d0
       IF( .NOT. lclassical) WRITE(iw,'(A,T71,F10.5)') &
           ' Maximum rotational displacement [degrees] = ',&
                           rmrot/pi*180.0d0
       WRITE(iw,*) '-----------------------------------------------',&
      '--------------------------------'
    ENDIF

    IF (moves%Quickstep%attempts .GT. 0) THEN
       WRITE(iw,'(A,T43,A)') " Quickstep Moves",&
          "Attempted       Accepted       Percent"
       WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%Quickstep%attempts,&
             moves%Quickstep%successes,&
             REAL(moves%Quickstep%successes,dp)/ &
             REAL(moves%Quickstep%attempts,dp)*100.0d0
       WRITE(iw,*) '-----------------------------------------------',&
      '--------------------------------'
    ENDIF

    IF (moves%swap%attempts .GT. 0) THEN
       WRITE(iw,'(A,T43,A)') " Swap Moves into this box",&
          "Attempted       Empty          Percent"
       WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%swap%attempts,&
             moves%empty,&
             REAL(moves%empty,dp)/ &
             REAL(moves%swap%attempts,dp)*100.0d0
       WRITE(iw,'(A,T43,A)') "                  Growths",&
          "Attempted       Sucessful      Percent"
       WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%swap%attempts,&
             moves%grown,&
             REAL(moves%grown,dp)/ &
             REAL(moves%swap%attempts,dp)*100.0d0
       WRITE(iw,'(A,T43,A)') "                    Total",&
          "Attempted       Accepted       Percent"
       WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%swap%attempts,&
             moves%swap%successes,&
             REAL(moves%swap%successes,dp)/ &
             REAL(moves%swap%attempts,dp)*100.0d0
       WRITE(iw,*) '-----------------------------------------------',&
      '--------------------------------'
    ENDIF

! now we write out information on the classical moves, if it's
! a classical simulations
    IF (lclassical) THEN
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
            REAL(move_updates%classical_bond%successes,dp)/ &
            REAL(move_updates%classical_bond%attempts,dp)*100.0d0
          WRITE(iw,'(A,T71,F10.5)') &
            ' Maximum bond displacement [angstroms]= ',&
            rmbond*angstrom
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
            REAL(move_updates%classical_angle%successes,dp)/ &
            REAL(move_updates%classical_angle%attempts,dp)*100.0d0
          WRITE(iw,'(A,T71,F10.5)') &
            ' Maximum angle displacement [degrees]= ',&
            rmangle/pi*180.0d0
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
            REAL(move_updates%classical_trans%successes,dp)/ &
            REAL(move_updates%classical_trans%attempts,dp)*100.0d0
          WRITE(iw,'(A,T71,F10.5)') &
         ' Maximum molecular translational displacement [angstroms]= ',&
                  rmtrans*angstrom
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
            REAL(move_updates%classical_rot%successes,dp)/ &
            REAL(move_updates%classical_rot%attempts,dp)*100.0d0
          WRITE(iw,'(A,T71,F10.5)') &
           ' Maximum rotational displacement [degrees] = ',&
                           rmrot/pi*180.0d0
          WRITE(iw,*) '-------------------------------------------',&
                   '------------------------------------'
       ENDIF

! now volume moves
       IF (moves%volume%attempts .GT. 0 .AND. nvolmoves .GT. 1)&
             THEN
          moves%classical_volume%attempts=&
             moves%classical_volume%attempts
          moves%classical_volume%successes=&
             moves%classical_volume%successes
          WRITE(iw,'(A,T43,A)') " Volume Changes",&
            "Attempted       Accepted       Percent"
          WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
            moves%classical_volume%attempts,&
            moves%classical_volume%successes,&
            REAL(moves%classical_volume%successes,dp)/ &
            REAL(moves%classical_volume%attempts,dp)*100.0d0
          WRITE(iw,'(A,T71,F10.5)') &
           ' Maximum Volume displacement [angstroms^3] = ',&
                     rmvolume*(angstrom)**3
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
!!****s* mc_misc.F/mc_make_dat_file
!!
!!   NAME
!!     mc_make_dat_file
!!
!!   FUNCTION
!!     writes a new input file that CP2K can read in for when we want
!!     to change a force env (change molecules or cell length)
!!
!!   ARGUMENTS
!!     - natoms: the total number of atoms
!!     - coordinates: the coordiantes of the atoms in the force_env (a.u.)
!!     - box_length: the length of the simulation box (angstrom)
!!     - filename: the name of the file to write to
!!     - force_env: the force environment that holds the parameters needed
!!                  to write on the file
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE mc_make_dat_file(coordinates,natoms,box_length,filename,force_env,&
                            fft_lib)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: coordinates
    INTEGER, INTENT(IN)                      :: natoms
    REAL(KIND=dp), INTENT(IN)                :: box_length
    CHARACTER(LEN=*), INTENT(IN)             :: filename
    CHARACTER(LEN=*), INTENT(IN)             :: fft_lib
    TYPE(force_env_type),POINTER             :: force_env

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_make_dat_file"

    TYPE(cell_type),POINTER                  :: cell_ref
    TYPE(qs_environment_type),POINTER        :: qs_env
    TYPE(dft_control_type),POINTER           :: dft_control
  
    INTEGER                                  :: handle, iatom, unit
    LOGICAL                                  :: use_ref_cell
    REAL(KIND=dp),DIMENSION(1:3)             :: abc
!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! grab some of the structures we're gonna need
      CALL force_env_get(force_env,qs_env=qs_env)
      CALL get_qs_env(qs_env,use_ref_cell=use_ref_cell,cell_ref=cell_ref,&
         dft_control=dft_control)
      CALL get_cell(cell_ref,abc=abc)

! open the file
      CALL open_file(file_name=filename,unit_number=unit,&
         file_action='WRITE',file_status='REPLACE')

! write all the stuff down
      WRITE(unit,10000) "&CP2K"
      WRITE(unit,10000) ' PROGRAM      Quickstep'
      WRITE(unit,10000) ' PROJECT      H2O_MC'
      WRITE(unit,10000) ' IOLEVEL      0'
      WRITE(unit,20001) ' FFTLIB       ',fft_lib
      WRITE(unit,10000) ' RUN_TYPE     MC'
      WRITE(unit,10000) "&END"
      WRITE(unit,*) 
      WRITE(unit,10000) "&PRINT low"
      WRITE(unit,10000) "&END"
      WRITE(unit,*) 
      WRITE(unit,10000) "&MC"
      WRITE(unit,10000) "&END"
      WRITE(unit,*) 
      WRITE(unit,10000) "&DFT"
!!!!!!!!!!!!!!!!!!!!!!
!      WRITE(unit,*) ' FUNCTIONAL   Pade'
      WRITE(unit,10000) ' FUNCTIONAL   BLYP'

!      WRITE(unit,*) ' XC_SMOOTH_RHO   ',&
!         ADJUSTR(xc_get_rho_smooth_label(dft_control%xc_rho_smooth_id))
      WRITE(unit,10000) ' XC_SMOOTH_RHO   NONE'
!      WRITE(unit,10000) ' XC_SMOOTH_RHO   NN50'
!      WRITE(unit,10000) ' XC_SMOOTH_RHO   NN10'
!      WRITE(unit,10000) ' XC_SMOOTH_RHO    SPLINE2'
!      WRITE(unit,10000) ' XC_SMOOTH_RHO    SPLINE3'
!     WRITE(unit,*) ' XC_DERIV     ',&
!         ADJUSTR(xc_get_deriv_method_label(dft_control%xc_deriv_method_id))
      WRITE(unit,10000) ' XC_DERIV      SPLINE2'
!      WRITE(unit,10000) ' XC_DERIV      SPLINE2_smooth'
!      WRITE(unit,10000) ' XC_DERIV      SPLINE3_smooth'
!      WRITE(unit,10000) ' XC_DERIV      NN10_smooth'
!      WRITE(unit,10000) ' XC_DERIV      NN50_smooth'
      WRITE(unit,10000) "&END"
      WRITE(unit,*) 
      WRITE(unit,10000) "&QS"

! I need to undo the unit change that's done in 
!       cp_control_units.F/read_qs_control
      WRITE(unit,20000) ' CUTOFF       ',&
         dft_control % qs_control % cutoff*2.0_dp
      WRITE(unit,10000) ' EXTRAPOLATION       USE_PREV_WF'
      WRITE(unit,10000) "&END"
      WRITE(unit,*) 
      WRITE(unit,10000) "&SCF"
      WRITE(unit,10000) ' GUESS        ATOMIC'
!      WRITE(unit,10000) ' GUESS        RESTART'
!      WRITE(unit,10000) ' OT'
      WRITE(unit,10000) "&END"
      WRITE(unit,*) 
!!!!!!!!!!!!!!!!!!!
!      WRITE(unit,10000) "&OTSCF"
!      WRITE(unit,10000) ' MINIMIZER DIIS 7'
!      WRITE(unit,10000) "&END"
!      WRITE(unit,*) 
!!!!!!!!!!!!!!!!!!!!!!!!!!
      WRITE(unit,10000) "&TOPOLOGY"
      WRITE(unit,10000) ' CONNECTIVITY   OHH_MC'
!      WRITE(unit,10000) ' BONDPARM COVALENT 1.6'
      WRITE(unit,10000) "&END"
      WRITE(unit,*) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      WRITE(unit,10000) "&KIND H"
!      WRITE(unit,10000) ' BASIS_SET    DZVP-GTH-PADE'
!      WRITE(unit,10000) ' POTENTIAL    GTH-PADE-q1'
      WRITE(unit,10000) ' BASIS_SET    TZV2P-GTH-BLYP'
      WRITE(unit,10000) ' POTENTIAL    GTH-BLYP-q1'
      WRITE(unit,10000) "&END"
      WRITE(unit,*)  
      WRITE(unit,10000) "&KIND O"
!      WRITE(unit,10000) ' BASIS_SET    DZVP-GTH-PADE'
!      WRITE(unit,10000) ' POTENTIAL    GTH-PADE-q6'
      WRITE(unit,10000) ' BASIS_SET    TZV2P-GTH-BLYP'
      WRITE(unit,10000) ' POTENTIAL    GTH-BLYP-q6'
      WRITE(unit,10000) "&END"
      WRITE(unit,*) 
      WRITE(unit,10000) "&CELL"
      WRITE(unit,10000) ' UNIT         ANGSTROM'
      WRITE(unit,40000) '  ABC          ',&
         box_length,box_length,box_length
      WRITE(unit,10000) "&END"
      WRITE(unit,*) 

! let's try to deal with the reference cell
      IF(use_ref_cell) THEN
         WRITE(unit,10000) '&CELL_REF'
         WRITE(unit,10000) ' UNIT         ANGSTROM'
         WRITE(unit,40000) '  ABC          ',&
            abc(1:3)*angstrom
         WRITE(unit,10000) "&END"
         WRITE(unit,*) 
      ENDIF

      WRITE(unit,10000) "&COORD"
      DO iatom=1,natoms,3
         WRITE(unit,50000) ' O  ',&
            coordinates(1:3,iatom)*angstrom,'  H2O'
         WRITE(unit,50000) ' H  ',&
            coordinates(1:3,iatom+1)*angstrom,'  H2O'
         WRITE(unit,50000) ' H  ',&
            coordinates(1:3,iatom+2)*angstrom,'  H2O'
      ENDDO
      WRITE(unit,10000) "&END"

! close the file
      CALL close_file(unit_number=unit)

! specify some formats
10000 FORMAT(A)
20000 FORMAT(A,F8.2)
20001 FORMAT(A,A)
40000 FORMAT(A,3(F12.8,2X))
50000 FORMAT(A,3(F18.12,2X),A)

! end the timing
      CALL timestop(0.0d0,handle)

!STOP

END SUBROUTINE mc_make_dat_file

END MODULE mc_misc

!******************************************************************************

