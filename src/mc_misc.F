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
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
!  USE cp_control_types,                ONLY: dft_control_type
!  USE cp_control_utils,                ONLY: xc_get_rho_smooth_label,&
!                                             xc_get_deriv_method_label
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE external_potential_types,        ONLY: elp_potential_type,&
                                             get_potential,&
                                             gth_potential_type
  USE fist_environment_types,          ONLY: fist_environment_type,&
                                             get_fist_env
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                             fist_nonbond_env_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type,&
                                             use_fist_force,&
                                             use_kg_force,&
                                             use_qs_force
  USE global_types,                    ONLY: SILENT,&
                                             LOW,&
                                             MEDIUM,&
                                             HIGH,&
                                             DEBUG
  USE input_constants,                 ONLY: do_ewald_ewald,&
                                             do_ewald_none,&
                                             do_ewald_pme,&
                                             do_ewald_spme
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_get_subs_vals2,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: pi
  USE mc_types,                        ONLY: get_mc_par,&
                                             mc_averages_type,&
                                             mc_moves_type,&
                                             mc_simulation_parameters_type
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE molecule_kind_types,             ONLY: atom_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type,&
                                             bond_type,&
                                             bend_type,&
                                             torsion_type
  USE pair_potential_types,            ONLY: pair_potential_pp_type
  USE physcon,                         ONLY: angstrom,&
                                             kelvin
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PUBLIC :: final_mc_write,mc_averages_create,mc_averages_release,&
            mc_make_dat_file

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mc_misc'

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
!!   Suitable for parallel.
!!
!!   AUTHOR
!!     MJM
!!
!!****************************************************************************

SUBROUTINE mc_averages_create ( averages  )


    TYPE(mc_averages_type), POINTER          :: averages

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mc_averages_create'

    INTEGER                                  :: handle, stat

!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! allocate all the structures...not sure why, but it won't work otherwise
      ALLOCATE (averages,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "averages",0)

      averages%ave_energy=0.0E0_dp
      averages%ave_energy_squared=0.0E0_dp
      averages%ave_volume=0.0E0_dp
      averages%bond1=0.0E0_dp
      averages%bond2=0.0E0_dp
      averages%angle=0.0E0_dp
      averages%molecules=0.0E0_dp

! end the timing
       CALL timestop(0.0E0_dp,handle)


END SUBROUTINE mc_averages_create

!******************************************************************************
!!****** mc_misc.F/mc_averages_release
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
!!   Suitable for parallel.
!!
!!   AUTHOR
!!     MJM
!!
!!****************************************************************************

SUBROUTINE mc_averages_release ( averages  )


    TYPE(mc_averages_type), POINTER          :: averages

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mc_averages_release'

    INTEGER                                  :: handle, stat

!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! deallocate 
      DEALLOCATE (averages,stat=stat)
      IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "averages")

      NULLIFY(averages)

! end the timing
       CALL timestop(0.0E0_dp,handle)


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
!!     - iw: the unit to write to
!!     - energy_check: the sum of the energy changes of each move
!!     - initial_energy: the initial unbiased energy of the system
!!     - final_energy: the final unbiased energy of the system
!!     - averages: the structure that holds computed average properites for
!!              the simulation
!!
!!   Only use in serial.
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

SUBROUTINE final_mc_write (mc_par,moves,iw,energy_check,initial_energy,&
                           final_energy,averages)


    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(mc_moves_type), POINTER             :: moves
    INTEGER, INTENT(IN)                      :: iw
    REAL(KIND=dp), INTENT(IN)                :: energy_check, initial_energy, &
                                                final_energy
    TYPE(mc_averages_type), POINTER          :: averages

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'final_mc_write'

    INTEGER                                  :: handle, nvolmoves
    LOGICAL                                  :: lbias
    REAL(KIND=dp)                            :: pmswap, rmangle, rmbond, &
                                                rmdihedral, rmrot, rmtrans, &
                                                rmvolume

!------------------------------------------------------------------------------
! begin the timing of the subroutine

    CALL timeset(routineN,'I','',handle)

    CALL get_mc_par(mc_par,pmswap=pmswap,rmvolume=rmvolume,&
      lbias=lbias,rmbond=rmbond,rmangle=rmangle,rmdihedral=rmdihedral,&
      rmtrans=rmtrans,rmrot=rmrot,nvolmoves=nvolmoves)

! write out some data averaged over the whole simulation
    WRITE(iw,*)
    WRITE(iw,*) '****************************************************',&
                   '***************************'
    WRITE(iw,'(A,T66,F15.8)') " Average Energy [Hartrees]:", &
                averages%ave_energy
!    WRITE(iw,'(A,T66,F15.8)') " Average O-H distance 1 [angstroms]:",&
!                   averages%bond1*angstrom
!    WRITE(iw,'(A,T66,F15.8)') " Average O-H distance 2 [angstroms]:",&
!                   averages%bond2*angstrom
!    WRITE(iw,'(A,T66,F15.8)') " Average H-O-H angle [degrees]:",&
!                   averages%angle
    IF (pmswap .GT. 0.0E0_dp) THEN
       WRITE(iw,'(A,T66,F15.8)') " Average number of molecules:",&
                   averages%molecules
    ENDIF
!    WRITE(iw,'(A,A,T68,F13.10)') " Average Density ",&
!      "[g/cm^3]:",1.0E0_dp/averages%ave_volume*averages%molecules/&
!      (angstrom)**3*(1.0e8_dp)**3/n_avogadro*18.0184E0_dp
    WRITE(iw,'(A,A,T65,F16.6)') " Average Volume ",&
               "[angstroms**3]:",averages%ave_volume*angstrom**3

    WRITE(iw,*)

! write out acceptance rates for the moves
    WRITE(iw,*) "Move Data"
    WRITE(iw,*) '----------------------------------------------------',&
                   '---------------------------'
! AVBMC moves
    IF (moves%avbmc_inin%attempts .GT. 0) THEN
       WRITE(iw,'(A,T43,A)') " AVBMC moves from in to in",&
          "Attempted       Accepted       Percent"
       WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%avbmc_inin%attempts,&
             moves%avbmc_inin%successes,&
             REAL(moves%avbmc_inin%successes,dp)/ &
             REAL(moves%avbmc_inin%attempts,dp)*100.0E0_dp
       WRITE(iw,*) '-----------------------------------------------',&
      '--------------------------------'
    ENDIF
    IF (moves%avbmc_inout%attempts .GT. 0) THEN
       WRITE(iw,'(A,T43,A)') " AVBMC moves from in to out",&
          "Attempted       Accepted       Percent"
       WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%avbmc_inout%attempts,&
             moves%avbmc_inout%successes,&
             REAL(moves%avbmc_inout%successes,dp)/ &
             REAL(moves%avbmc_inout%attempts,dp)*100.0E0_dp
       WRITE(iw,*) '-----------------------------------------------',&
      '--------------------------------'
    ENDIF
    IF (moves%avbmc_outin%attempts .GT. 0) THEN
       WRITE(iw,'(A,T43,A)') " AVBMC moves from out to in",&
          "Attempted       Accepted       Percent"
       WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%avbmc_outin%attempts,&
             moves%avbmc_outin%successes,&
             REAL(moves%avbmc_outin%successes,dp)/ &
             REAL(moves%avbmc_outin%attempts,dp)*100.0E0_dp
       WRITE(iw,*) '-----------------------------------------------',&
      '--------------------------------'
    ENDIF
    IF (moves%avbmc_outout%attempts .GT. 0) THEN
       WRITE(iw,'(A,T43,A)') " AVBMC moves from out to out",&
          "Attempted       Accepted       Percent"
       WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%avbmc_outout%attempts,&
             moves%avbmc_outout%successes,&
             REAL(moves%avbmc_outout%successes,dp)/ &
             REAL(moves%avbmc_outout%attempts,dp)*100.0E0_dp
       WRITE(iw,*) '-----------------------------------------------',&
      '--------------------------------'
    ENDIF

! volume moves
    IF (moves%volume%attempts .GT. 0) THEN
       WRITE(iw,'(A,T43,A)') " Volume Moves",&
          "Attempted       Accepted       Percent"
       WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%volume%attempts,&
             moves%volume%successes,&
             REAL(moves%volume%successes,dp)/ &
             REAL(moves%volume%attempts,dp)*100.0E0_dp
       WRITE(iw,'(A,T71,F10.5)') &
           ' Maximum Volume displacement [angstroms^3] = ',&
                     rmvolume*(angstrom)**3
       WRITE(iw,*) '-----------------------------------------------',&
      '--------------------------------'
    ENDIF

    IF (moves%angle%attempts .GT. 0 .OR. &
      moves%bond%attempts .GT. 0 .OR. moves%dihedral%attempts .GT. 0 ) THEN
       WRITE(iw,'(A,T43,A)') " Conformational Moves",&
          "Attempted       Accepted       Percent"
       WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
          moves%bond%attempts+moves%angle%attempts+moves%dihedral%attempts,&
          moves%bond%successes+moves%angle%successes+moves%dihedral%successes,&
          REAL(moves%bond%successes+moves%angle%successes+moves%dihedral%successes,dp)/ &
          REAL(moves%bond%attempts+moves%angle%attempts+moves%dihedral%attempts,dp)*100.0E0_dp
       IF (moves%bond%attempts .GT. 0 ) THEN
          WRITE(iw,'(A,T43,A)') " Bond Changes",&
            "Attempted       Accepted       Percent"
          WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
            moves%bond%attempts,&
            moves%bond%successes,&
            REAL(moves%bond%successes,dp)/ &
            REAL(moves%bond%attempts,dp)*100.0E0_dp
          IF( .NOT. lbias) WRITE(iw,'(A,T71,F10.5)') &
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
            REAL(moves%angle%attempts,dp)*100.0E0_dp
          IF( .NOT. lbias) WRITE(iw,'(A,T71,F10.5)') &
            ' Maximum angle displacement [degrees]= ',&
            rmangle/pi*180.0E0_dp
       ENDIF
       IF (moves%dihedral%attempts .GT. 0 ) THEN
          WRITE(iw,'(A,T43,A)') " Dihedral Changes",&
            "Attempted       Accepted       Percent"
          WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
            moves%dihedral%attempts,&
            moves%dihedral%successes,&
            REAL(moves%dihedral%successes,dp)/ &
            REAL(moves%dihedral%attempts,dp)*100.0E0_dp
          IF( .NOT. lbias) WRITE(iw,'(A,T71,F10.5)') &
            ' Maximum dihedral displacement [degrees]= ',&
            rmdihedral/pi*180.0E0_dp
       ENDIF
       WRITE(iw,'(A,A,I5)') " Conformational Moves Rejected Because",&
          "Box Was Empty: ",moves%empty_conf
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
             REAL(moves%trans%attempts,dp)*100.0E0_dp
       IF( .NOT. lbias) WRITE(iw,'(A,T71,F10.5)') &
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
             REAL(moves%rot%attempts,dp)*100.0E0_dp
       IF( .NOT. lbias) WRITE(iw,'(A,T71,F10.5)') &
           ' Maximum rotational displacement [degrees] = ',&
                           rmrot/pi*180.0E0_dp
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
             REAL(moves%Quickstep%attempts,dp)*100.0E0_dp
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
             REAL(moves%swap%attempts,dp)*100.0E0_dp
       WRITE(iw,'(A,T43,A)') "                  Growths",&
          "Attempted       Sucessful      Percent"
       WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%swap%attempts,&
             moves%grown,&
             REAL(moves%grown,dp)/ &
             REAL(moves%swap%attempts,dp)*100.0E0_dp
       WRITE(iw,'(A,T43,A)') "                    Total",&
          "Attempted       Accepted       Percent"
       WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
             moves%swap%attempts,&
             moves%swap%successes,&
             REAL(moves%swap%successes,dp)/ &
             REAL(moves%swap%attempts,dp)*100.0E0_dp
       WRITE(iw,*) '-----------------------------------------------',&
      '--------------------------------'
    ENDIF

! now we write out information on the classical moves, if it's
! a classical simulations
    IF (lbias) THEN
       WRITE(iw,*) "Biased Move Data"
       WRITE(iw,*) '------------------------------------------------',&
                   '-------------------------------'
! first, bond changes...record all the moves done since the last update,
! and print out the totals
       IF (moves%bond%attempts .GT. 0 ) THEN
          WRITE(iw,'(A,T43,A)') " Bond Changes",&
            "Attempted       Accepted       Percent"
          WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
            moves%bias_bond%attempts,&
            moves%bias_bond%successes,&
            REAL(moves%bias_bond%successes,dp)/ &
            REAL(moves%bias_bond%attempts,dp)*100.0E0_dp
          WRITE(iw,'(A,T71,F10.5)') &
            ' Maximum bond displacement [angstroms]= ',&
            rmbond*angstrom
          WRITE(iw,*) '-------------------------------------------',&
                   '------------------------------------'
       ENDIF

! now angle changes
       IF (moves%angle%attempts .GT. 0 ) THEN
          WRITE(iw,'(A,T43,A)') " Angle Changes",&
            "Attempted       Accepted       Percent"
          WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
            moves%bias_angle%attempts,&
            moves%bias_angle%successes,&
            REAL(moves%bias_angle%successes,dp)/ &
            REAL(moves%bias_angle%attempts,dp)*100.0E0_dp
          WRITE(iw,'(A,T71,F10.5)') &
            ' Maximum angle displacement [degrees]= ',&
            rmangle/pi*180.0E0_dp
          WRITE(iw,*) '-------------------------------------------',&
                   '------------------------------------'
       ENDIF

! now dihedral changes
       IF (moves%dihedral%attempts .GT. 0 ) THEN
          WRITE(iw,'(A,T43,A)') " Dihedral Changes",&
            "Attempted       Accepted       Percent"
          WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
            moves%bias_dihedral%attempts,&
            moves%bias_dihedral%successes,&
            REAL(moves%bias_dihedral%successes,dp)/ &
            REAL(moves%bias_dihedral%attempts,dp)*100.0E0_dp
          WRITE(iw,'(A,T71,F10.5)') &
            ' Maximum dihedral displacement [degrees]= ',&
            rmdihedral/pi*180.0E0_dp
          WRITE(iw,*) '-------------------------------------------',&
                   '------------------------------------'
       ENDIF

! now molecule translations
       IF (moves%trans%attempts .GT. 0 ) THEN
          WRITE(iw,'(A,T43,A)') " Molecule Translations",&
            "Attempted       Accepted       Percent"
          WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
            moves%bias_trans%attempts,&
            moves%bias_trans%successes,&
            REAL(moves%bias_trans%successes,dp)/ &
            REAL(moves%bias_trans%attempts,dp)*100.0E0_dp
          WRITE(iw,'(A,T71,F10.5)') &
         ' Maximum molecular translational displacement [angstroms]= ',&
                  rmtrans*angstrom
          WRITE(iw,*) '-------------------------------------------',&
                   '------------------------------------'
       ENDIF

! now molecule rotations
       IF (moves%rot%attempts .GT. 0 ) THEN
          WRITE(iw,'(A,T43,A)') " Molecule Rotations",&
            "Attempted       Accepted       Percent"
          WRITE(iw,'(T46,I6,9X,I6,7X,F7.3)') &
            moves%bias_rot%attempts,&
            moves%bias_rot%successes,&
            REAL(moves%bias_rot%successes,dp)/ &
            REAL(moves%bias_rot%attempts,dp)*100.0E0_dp
          WRITE(iw,'(A,T71,F10.5)') &
           ' Maximum rotational displacement [degrees] = ',&
                           rmrot/pi*180.0E0_dp
          WRITE(iw,*) '-------------------------------------------',&
                   '------------------------------------'
       ENDIF

    ENDIF



! see if the energies add up properly
    IF(ABS(initial_energy+energy_check-final_energy) .GT. 0.0000001E0_dp) &
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
    CALL timestop(0.0E0_dp,handle)

 
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
!!     - coordinates: the coordiantes of the atoms in the force_env (a.u.)
!!     - natoms_tot: the total number of atoms
!!     - box_length: the length of all sides of the simulation box (angstrom)
!!     - filename: the name of the file to write to
!!     - force_env: the force environment that holds the parameters needed
!!                  to write on the file
!!     - fft_lib: the FFT library we're using (FFTW,FFTSG)
!!
!!   Only use in serial.
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE mc_make_dat_file(coordinates,natoms_tot,box_length,filename,&
      force_env,fft_lib,print_level)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: coordinates
    INTEGER, INTENT(IN)                      :: natoms_tot,print_level
    REAL(KIND=dp), DIMENSION(1:3), &
      INTENT(IN)                             :: box_length
    CHARACTER(LEN=*), INTENT(IN)             :: filename
    TYPE(force_env_type), POINTER            :: force_env
    CHARACTER(LEN=*), INTENT(IN)             :: fft_lib

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mc_make_dat_file', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=5)                         :: cfile
    CHARACTER(LEN=default_string_length)     :: c_val, name, unit_type,&
         name_a,name_b,name_c,name_d
    INTEGER :: ewald_type, gmax, handle, i_val, iatom, ichain, ifunct, itype, &
      jtype, natoms, nchains, ns_max, ntypes, o_spline, unit,&
      nbond,nbend,ntorsion,ibond,ibend,itorsion,jbond,stat,current_bond,&
      jbend,current_bend,imul,current_torsion,jtorsion
    LOGICAL                                  :: failure, use_ref_cell
    REAL(KIND=dp)                            :: alpha, charge, epsilon, r_val
    REAL(KIND=dp), DIMENSION(1:3)            :: abc
    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(cp_error_type)                      :: error
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(elp_potential_type), POINTER        :: elp_potential
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(fist_environment_type), POINTER     :: fist_env
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: input, xc_fun, xc_fun_section
    TYPE(bond_type), DIMENSION(:), POINTER   :: bond_list
    TYPE(bend_type), DIMENSION(:), POINTER   :: bend_list
    TYPE(torsion_type), DIMENSION(:), POINTER   :: torsion_list
!    TYPE(input_info_type), POINTER           :: inp_info
    CHARACTER(default_string_length),DIMENSION(:),ALLOCATABLE :: atom_a_done,&
         atom_b_done,atom_c_done,atom_d_done
!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! nullify some stuff
      NULLIFY(qs_env,atomic_kind_set,atomic_kind,gth_potential,&
         orb_basis_set,elp_potential,input,xc_fun,kg_env,&
         bond_list,bend_list,torsion_list)

! grab some of the structures we're gonna need
      SELECT CASE ( force_env%in_use )
      CASE ( use_fist_force )
         CALL force_env_get(force_env,fist_env=fist_env,subsys=subsys)
         CALL get_fist_env(fist_env,cell=cell,ewald_env=ewald_env,&
            atomic_kind_set=atomic_kind_set,&
            fist_nonbond_env=fist_nonbond_env,error=error)
         CALL ewald_env_get(ewald_env,unit_type=unit_type,&
            ewald_type=ewald_type,alpha=alpha,gmax=gmax,o_spline=o_spline,&
            ns_max=ns_max,epsilon=epsilon)
         CALL get_cell(cell,abc=abc)
      CASE ( use_qs_force )
         CALL force_env_get(force_env,qs_env=qs_env,subsys=subsys)
         CALL get_qs_env(qs_env,use_ref_cell=use_ref_cell,cell_ref=cell_ref,&
           input=input,atomic_kind_set=atomic_kind_set)
         CALL get_cell(cell_ref,abc=abc)

      CASE ( use_kg_force )
         CALL force_env_get(force_env,kg_env=kg_env,subsys=subsys)
         CALL get_kg_env(kg_env,use_ref_cell=use_ref_cell,cell_ref=cell_ref,&
           input=input,atomic_kind_set=atomic_kind_set)
         CALL get_cell(cell_ref,abc=abc)

      CASE default
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
      END SELECT

! get some molecule information
      CALL cp_subsys_get(subsys(1)%subsys, &
         molecule_kinds_new=molecule_kinds_new,error=error)
      molecule_kind => molecule_kinds_new%els(1)
      CALL get_molecule_kind(molecule_kind,atom_list=atom_list,&
         natom=natoms,bond_list=bond_list,bend_list=bend_list,&
         torsion_list=torsion_list,nbond=nbond,nbend=nbend,ntorsion=ntorsion)
      
! find out how many atom types we have
      ntypes=SIZE(atomic_kind_set(:))

! a way to check and see what bond, bend, and tors types we've
! already written
!      ALLOCATE (lbond(1:nbonmolecule_kind_set(nmolecule_kind),STAT=istat)
!    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
!                                     "molecule_kind_set",&
!                                     nmolecule_kind*int_size)

! open the file
      CALL open_file(file_name=filename,unit_number=unit,&
         file_action='WRITE',file_status='REPLACE')

! write all the stuff down
      WRITE(unit,10000) "&FORCE_EVAL"
      SELECT CASE ( force_env%in_use )
      CASE ( use_fist_force )
! grab some of the Ewald keywords
         WRITE(unit,10000) "  &MM"
         WRITE(unit,10000) "    &POISSON"
         WRITE(unit,10000) "      &EWALD"
         CALL convert_to_cp2k_units ( 'BOHR', length = alpha,&
            l_power = -1 )  
         SELECT CASE ( ewald_type )
         CASE ( do_ewald_ewald )
            WRITE(unit,10000) "       EWALD_TYPE ewald"
            WRITE(unit,20000) "       ALPHA ",alpha
            WRITE(unit,20003) "       GMAX ",gmax
         CASE ( do_ewald_none )
            WRITE(unit,10000) "       EWALD_TYPE none"
            WRITE(unit,20000) "       ALPHA ",alpha
            WRITE(unit,20003) "       GMAX ",gmax
         CASE ( do_ewald_pme )
            WRITE(unit,10000) "       EWALD_TYPE pme"
            WRITE(unit,20000) "       ALPHA ",alpha
            WRITE(unit,20003) "       NS_MAX ",ns_max
            WRITE(unit,20000) "       EPSILON ",epsilon
         CASE ( do_ewald_spme )
            WRITE(unit,10000) "       EWALD_TYPE pme"
            WRITE(unit,20000) "       ALPHA ",alpha
            WRITE(unit,20003) "       GMAX ",gmax
            WRITE(unit,20003) "       O_SPLINE ",o_spline
         END SELECT
         WRITE(unit,10000) "      &END EWALD"
         WRITE(unit,10000) "    &END POISSON"


! rcut
         CALL fist_nonbond_env_get (fist_nonbond_env, potparm=potparm)

! now some of the force field stuff
         WRITE(unit,10000) "    &FORCEFIELD"
! doesn't look like these values are stored past read-in
!         WRITE(unit,10000) "      parmfile water.pot"
!         WRITE(unit,10000) "      parmtype CHM"

         
! charges
         DO itype=1,ntypes
            atomic_kind => atomic_kind_set(itype)
            CALL get_atomic_kind(atomic_kind=atomic_kind,&
              name=name,elp_potential=elp_potential)
            CALL get_potential(potential=elp_potential,&
              qeff=charge)
            WRITE(unit,10000) "      &CHARGE"
            WRITE(unit,20001) "        ATOM ",TRIM(ADJUSTL(name))
            WRITE(unit,20002) "        CHARGE ",charge
            WRITE(unit,10000) "      &END CHARGE"
         ENDDO

! need to find all the bonds and write them out
         ALLOCATE(atom_a_done(1:nbond),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(atom_b_done(1:nbond),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         atom_a_done(:)=""
         atom_b_done(:)=""
         current_bond=1
         outer_bond:DO ibond=1,nbond
            name_a=atom_list(bond_list(ibond)%a)%name
            name_b=atom_list(bond_list(ibond)%b)%name
            DO jbond=1,current_bond
               IF(atom_a_done(jbond)==name_a .AND. atom_b_done(jbond)==name_b&
                    .OR. atom_a_done(jbond)==name_b .AND. atom_b_done(jbond)==name_a)THEN
                  CYCLE outer_bond
               ENDIF
            ENDDO
            atom_a_done(current_bond)=name_a
            atom_b_done(current_bond)=name_b
            current_bond=current_bond+1
            WRITE(unit,10000) "      &BOND"
            WRITE(unit,30002) "        ATOMS ",TRIM(ADJUSTL(name_a)),TRIM(ADJUSTL(name_b))
            WRITE(unit,20002) "        K ",bond_list(ibond)%bond_kind%k(1)
            WRITE(unit,20002) "        R0 ",bond_list(ibond)%bond_kind%r0
            WRITE(unit,10000) "      &END BOND"
         ENDDO outer_bond
         DEALLOCATE(atom_a_done)
         DEALLOCATE(atom_b_done)

! same thing for bends
         ALLOCATE(atom_a_done(1:nbend),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(atom_b_done(1:nbend),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(atom_c_done(1:nbend),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         atom_a_done(:)=""
         atom_b_done(:)=""
         atom_c_done(:)=""
         current_bend=1
         outer_bend:DO ibend=1,nbend
            name_a=atom_list(bend_list(ibend)%a)%name
            name_b=atom_list(bend_list(ibend)%b)%name
            name_c=atom_list(bend_list(ibend)%c)%name
            DO jbend=1,current_bend
               IF(atom_a_done(jbend)==name_a .AND. atom_b_done(jbend)==name_b &
                    .AND.atom_c_done(jbend)==name_c &
                    .OR. atom_a_done(jbend)==name_c .AND. atom_b_done(jbend)==name_b .AND. &
                    atom_c_done(jbend)==name_a)THEN
                  CYCLE outer_bend
               ENDIF
            ENDDO
            atom_a_done(current_bend)=name_a
            atom_b_done(current_bend)=name_b
            atom_c_done(current_bend)=name_c
            current_bend=current_bend+1
            WRITE(unit,10000) "      &BEND"
            WRITE(unit,40004) "        ATOMS ",TRIM(ADJUSTL(name_a)),TRIM(ADJUSTL(name_b)),TRIM(ADJUSTL(name_c))
            WRITE(unit,20002) "        K ",bend_list(ibend)%bend_kind%k
            WRITE(unit,20002) "        THETA0 ",bend_list(ibend)%bend_kind%theta0
            WRITE(unit,10000) "      &END BEND"
         ENDDO outer_bend
         DEALLOCATE(atom_a_done)
         DEALLOCATE(atom_b_done)
         DEALLOCATE(atom_c_done)

! finally, write out the torsions
         ALLOCATE(atom_a_done(1:ntorsion),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(atom_b_done(1:ntorsion),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(atom_c_done(1:ntorsion),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(atom_d_done(1:ntorsion),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         atom_a_done(:)=""
         atom_b_done(:)=""
         atom_c_done(:)=""
         atom_d_done(:)=""
         current_torsion=1
         outer_torsion:DO itorsion=1,ntorsion
            name_a=atom_list(torsion_list(itorsion)%a)%name
            name_b=atom_list(torsion_list(itorsion)%b)%name
            name_c=atom_list(torsion_list(itorsion)%c)%name
            name_d=atom_list(torsion_list(itorsion)%d)%name
            DO jtorsion=1,current_torsion
               IF(atom_a_done(jtorsion)==name_a .AND. atom_b_done(jtorsion)==name_b &
                    .AND.atom_c_done(jtorsion)==name_c .AND. atom_d_done(jtorsion)==name_d &
                    .OR. atom_a_done(jtorsion)==name_d .AND. atom_b_done(jtorsion)==name_c .AND. &
                    atom_c_done(jtorsion)==name_b .AND. atom_d_done(jtorsion)==name_a)THEN
                  CYCLE outer_torsion
               ENDIF
            ENDDO
            atom_a_done(current_torsion)=name_a
            atom_b_done(current_torsion)=name_b
            atom_c_done(current_torsion)=name_c
            atom_d_done(current_torsion)=name_d
            current_torsion=current_torsion+1
            DO imul=1,torsion_list(itorsion)%torsion_kind%nmul
               WRITE(unit,10000) "      &TORSION"
               WRITE(unit,50001) "        ATOMS ",TRIM(ADJUSTL(name_a)),TRIM(ADJUSTL(name_b)),&
                    TRIM(ADJUSTL(name_c)),TRIM(ADJUSTL(name_d))
               WRITE(unit,20003) "        M ",torsion_list(itorsion)%torsion_kind%m(imul)
               WRITE(unit,20002) "        K ",torsion_list(itorsion)%torsion_kind%k(imul)
               WRITE(unit,20002) "        PHI0 ",torsion_list(itorsion)%torsion_kind%cosphi0(imul)
               WRITE(unit,10000) "      &END TORSION"
            ENDDO
               
         ENDDO outer_torsion
         DEALLOCATE(atom_a_done)
         DEALLOCATE(atom_b_done)
         DEALLOCATE(atom_c_done)
         DEALLOCATE(atom_d_done)


! nonbonded...only works for LJ at the moment
         WRITE(unit,10000) "      &NONBONDED"
         DO itype=1,SIZE(potparm%pot,1)
            DO jtype=itype,SIZE(potparm%pot,2)
               WRITE(unit,10000) "        &LENNARD-JONES"
               WRITE(unit,30002) "          ATOMS ",&
                  TRIM(ADJUSTL(potparm%pot(itype,jtype)%pot%at1)),&
                  TRIM(ADJUSTL(potparm%pot(itype,jtype)%pot%at2))
               WRITE(unit,20002) "          EPSILON ",&
                  potparm%pot(itype,jtype)%pot%lj%epsilon*kelvin
               WRITE(unit,20002) "          SIGMA ",&
                  potparm%pot(itype,jtype)%pot%lj%sigma6**(1.0_dp/6.0_dp)*&
                  angstrom
               WRITE(unit,20000) "          RCUT ",&
                  potparm%pot(itype,jtype)%pot%rcutsq**0.5_dp*angstrom
               WRITE(unit,10000) "        &END LENNARD-JONES"
            ENDDO
         ENDDO
         WRITE(unit,10000) "      &END NONBONDED"
         WRITE(unit,10000) "    &END FORCEFIELD"
         WRITE(unit,10000) "  &END MM"


      CASE ( use_qs_force )
         WRITE(unit,10000) "  &DFT"
         WRITE(unit,10000) "    &MGRID"
         CALL section_vals_val_get(input,"DFT%MGRID%CUTOFF",r_val=r_val)
         WRITE(unit,20000) "      CUTOFF  ",2.0E0_dp*r_val
         CALL section_vals_val_get(input,"DFT%MGRID%NGRIDS",i_val=i_val)
         WRITE(unit,20003) "      NGRIDS  ",i_val
         WRITE(unit,10000) "    &END MGRID"
         WRITE(unit,10000) "    &QS"
         WRITE(unit,10000) "      EXTRAPOLATION USE_PREV_WF"
         WRITE(unit,10000) "    &END QS"
         WRITE(unit,10000) "    &SCF"
         CALL section_vals_val_get(input,"DFT%SCF%SCF_GUESS",c_val=c_val)
         WRITE(unit,20001) "      SCF_GUESS  ",TRIM(ADJUSTL(c_val))
         WRITE(unit,10000) "    &END SCF"
         WRITE(unit,10000) "    &XC"
         xc_fun_section => section_vals_get_subs_vals(input,&
            "DFT%XC%XC_FUNCTIONAL")
         WRITE(unit,10000) "      &XC_FUNCTIONAL"
         DO ifunct=1,4
            xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifunct)
            IF(ASSOCIATED(xc_fun)) THEN
               WRITE(unit,20001) "        &",TRIM(xc_fun%section%name)
               WRITE(unit,20001) "        &END ",TRIM(xc_fun%section%name)
            ENDIF
         ENDDO
         WRITE(unit,10000) "      &END XC_FUNCTIONAL"
         WRITE(unit,10000) "      &XC_GRID"
         CALL section_vals_val_get(input,"DFT%XC%XC_GRID%XC_SMOOTH_RHO",&
           c_val=c_val)
         WRITE(unit,20001) "        XC_SMOOTH_RHO ",TRIM(ADJUSTL(c_val))
         CALL section_vals_val_get(input,"DFT%XC%XC_GRID%XC_DERIV",&
           c_val=c_val)
         WRITE(unit,20001) "        XC_DERIV ",TRIM(ADJUSTL(c_val))
         WRITE(unit,10000) "      &END XC_GRID"
         WRITE(unit,10000) "    &END XC"
         WRITE(unit,10000) "  &END DFT"
      CASE ( use_kg_force )
         WRITE(unit,10000) "  &DFT"
         WRITE(unit,10000) "    &MGRID"
         CALL section_vals_val_get(input,"DFT%MGRID%CUTOFF",r_val=r_val)
         WRITE(unit,20000) "      CUTOFF  ",2.0E0_dp*r_val
         CALL section_vals_val_get(input,"DFT%MGRID%NGRIDS",i_val=i_val)
         WRITE(unit,20003) "      NGRIDS  ",i_val
         WRITE(unit,10000) "    &END MGRID"
         WRITE(unit,10000) "    &QS"
         CALL section_vals_val_get(input,"DFT%QS%EXTRAPOLATION",c_val=c_val)
         WRITE(unit,20001) "      EXTRAPOLATION  ",TRIM(ADJUSTL(c_val))
         CALL section_vals_val_get(input,"DFT%QS%EPS_CORE_CHARGE",r_val=r_val)
         WRITE(unit,20004) "      EPS_CORE_CHARGE  ",r_val
         CALL section_vals_val_get(input,"DFT%QS%EPS_PPL",r_val=r_val)
         WRITE(unit,20004) "      EPS_PPL  ",r_val
         CALL section_vals_val_get(input,"DFT%QS%EPS_PGF_ORB",r_val=r_val)
         WRITE(unit,20004) "      EPS_PGF_ORB  ",r_val
         CALL section_vals_val_get(input,"DFT%QS%EPS_GVG_RSPACE",r_val=r_val)
         WRITE(unit,20004) "      EPS_GVG_RSPACE  ",r_val
         CALL section_vals_val_get(input,"DFT%QS%EPS_RHO",r_val=r_val)
         WRITE(unit,20004) "      EPS_RHO  ",r_val
         CALL section_vals_val_get(input,"DFT%QS%METHOD",c_val=c_val)
         WRITE(unit,20001) "      METHOD  ",TRIM(ADJUSTL(c_val))
         WRITE(unit,10000) "    &END QS"
         WRITE(unit,10000) "    &SCF"
         CALL section_vals_val_get(input,"DFT%SCF%SCF_GUESS",c_val=c_val)
         WRITE(unit,20001) "      SCF_GUESS  ",TRIM(ADJUSTL(c_val))
         WRITE(unit,10000) "    &END SCF"
         WRITE(unit,10000) "    &XC"
         xc_fun_section => section_vals_get_subs_vals(input,&
            "DFT%XC%XC_FUNCTIONAL")
! take care of the functionals...I had to put the extra if statement in for the kinetic
! energy functional, as that doesn't just use the section name
         WRITE(unit,10000) "      &XC_FUNCTIONAL"
         DO ifunct=1,4
            xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifunct)
            IF(ASSOCIATED(xc_fun)) THEN
               WRITE(unit,20001) "        &",TRIM(xc_fun%section%name)
               IF(TRIM(xc_fun%section%name) == "KE_GGA") THEN
                  CALL section_vals_val_get(input,"DFT%XC%XC_FUNCTIONAL%KE_GGA%FUNCTIONAL",&
                       c_val=c_val)
                  WRITE(unit,20001) "          FUNCTIONAL ",TRIM(c_val)
               ENDIF
               WRITE(unit,20001) "        &END ",TRIM(xc_fun%section%name)
            ENDIF
         ENDDO
         WRITE(unit,10000) "      &END XC_FUNCTIONAL"
         WRITE(unit,10000) "      &XC_GRID"
         CALL section_vals_val_get(input,"DFT%XC%XC_GRID%XC_SMOOTH_RHO",&
           c_val=c_val)
         WRITE(unit,20001) "        XC_SMOOTH_RHO ",TRIM(ADJUSTL(c_val))
         CALL section_vals_val_get(input,"DFT%XC%XC_GRID%XC_DERIV",&
           c_val=c_val)
         WRITE(unit,20001) "        XC_DERIV ",TRIM(ADJUSTL(c_val))
         WRITE(unit,10000) "      &END XC_GRID"
         WRITE(unit,10000) "    &END XC"
         WRITE(unit,10000) "  &END DFT"
      CASE default
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
      END SELECT

      WRITE(unit,10000) "  &SUBSYS"
! cell information
      WRITE(unit,10000) "    &CELL"
      WRITE(unit,40000) '      ABC  ',box_length(1:3)
      WRITE(unit,10000) '      UNIT ANGSTROM'
      IF(force_env%in_use == use_qs_force) THEN
         IF(use_ref_cell) THEN
            WRITE(unit,10000) '       &CELL_REF'
            WRITE(unit,40000) '         ABC ',&
               abc(1:3)*angstrom
            WRITE(unit,10000) '         UNIT ANGSTROM'
            WRITE(unit,10000) '       &END CELL_REF'
         ENDIF
      ENDIF
      WRITE(unit,10000) '    &END CELL'


! write out the coordinate information...currently only works for
! single component systems
      WRITE(unit,10000) '    &COORD'
      IF(MOD(natoms_tot,natoms) .NE. 0) THEN
         CALL stop_program("mc_make_dat_file",&
              "Not a single component system?")
      ENDIF
      nchains=natoms_tot/natoms

      IF(nchains .NE. 0) THEN
         DO ichain=1,nchains
            DO iatom=1,natoms
               WRITE(unit,40003) TRIM(ADJUSTL(atom_list(iatom)%name)),&
               coordinates(1:3,(ichain-1)*natoms+iatom)*angstrom
            ENDDO
         ENDDO
      ELSE
         DO iatom=1,natoms
            WRITE(unit,40003) TRIM(ADJUSTL(atom_list(iatom)%name)),&
               REAL(iatom,dp),REAL(iatom,dp),REAL(iatom,dp)
         ENDDO
      ENDIF
      WRITE(unit,10000) '    &END COORD'

! now the kind types, for the basis sets and potentials
      IF(force_env%in_use == use_qs_force .OR. &
         force_env%in_use == use_kg_force)THEN
         DO itype=1,ntypes
            atomic_kind => atomic_kind_set(itype)
            CALL get_atomic_kind(atomic_kind=atomic_kind,&
            element_symbol=element_symbol,gth_potential=gth_potential,&
            orb_basis_set=orb_basis_set,elp_potential=elp_potential)
            
            CALL get_gto_basis_set(orb_basis_set,name=name)
            
            WRITE(unit,20001) "    &KIND ",element_symbol
            WRITE(unit,20001) "      BASIS_SET    ",TRIM(ADJUSTL(name))
            
            IF(ASSOCIATED(gth_potential)) THEN
               CALL get_potential(gth_potential,name=name)
               WRITE(unit,20001) "      POTENTIAL    ",TRIM(ADJUSTL(name))
            ELSE
               CALL get_potential(elp_potential,name=name)
               WRITE(unit,20001) "      POTENTIAL    ",TRIM(ADJUSTL(name))
            ENDIF
            WRITE(unit,10000) "    &END KIND"
         ENDDO      
      ENDIF
! write the topology section
      WRITE(unit,10000) "    &TOPOLOGY"
      WRITE(unit,10000) "      CONNECTIVITY PSF"
      WRITE(cfile,'(I4)') natoms_tot/natoms
      IF(nchains == 0) cfile="1"
      SELECT CASE ( force_env%in_use )
      CASE ( use_fist_force )
         WRITE(unit,10000) '      CONN_FILE topology_fist_' &
              // TRIM(ADJUSTL(cfile)) // '.psf'
      CASE ( use_qs_force )
         WRITE(unit,10000) '      CONN_FILE topology_atoms_' &
              // TRIM(ADJUSTL(cfile)) // '.psf'
      CASE ( use_kg_force )
         WRITE(unit,10000) '      CONN_FILE topology_atoms_' &
              // TRIM(ADJUSTL(cfile)) // '.psf'
      CASE default
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
      END SELECT
      WRITE(unit,10000) "    &END TOPOLOGY"
      WRITE(unit,10000) "  &END SUBSYS"
      WRITE(unit,10000) "&END FORCE_EVAL"
      WRITE(unit,10000) "&GLOBAL"
      SELECT CASE ( force_env%in_use )
      CASE ( use_fist_force )
         WRITE(unit,10000) ' PROGRAM      FIST'
      CASE ( use_qs_force )
         WRITE(unit,10000) ' PROGRAM      Quickstep'
      CASE ( use_kg_force )
         WRITE(unit,10000) ' PROGRAM      KG'
      CASE default
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
      END SELECT
      WRITE(unit,10000) ' PROJECT      H2O_MC'
      WRITE(unit,20001) ' FFTLIB       ',TRIM(ADJUSTL(fft_lib))
! need to write the correct print level
      SELECT CASE (print_level)
      CASE(SILENT)
         WRITE(unit,10000) ' PRINT_LEVEL SILENT'
      CASE(LOW)
         WRITE(unit,10000) ' PRINT_LEVEL LOW'
      CASE(MEDIUM)
         WRITE(unit,10000) ' PRINT_LEVEL MEDIUM'
      CASE(HIGH)
         WRITE(unit,10000) ' PRINT_LEVEL HIGH'
      CASE(DEBUG)
         WRITE(unit,10000) ' PRINT_LEVEL DEBUG'
      END SELECT
      WRITE(unit,10000) "&END"

! close the file
      CALL close_file(unit_number=unit)

! specify some formats
10000 FORMAT(A)
20000 FORMAT(A,F8.2)
20001 FORMAT(A,A)
20002 FORMAT(A,F12.6)
20003 FORMAT(A,I10)
20004 FORMAT(A,E12.6)
30000 FORMAT(A,F12.8,2X,I8)
30001 FORMAT(A,A,F16.10)
30002 FORMAT(A,A,1X,A)
40000 FORMAT(A,3(F20.12,2X))
40001 FORMAT(A,F12.8,2X,I8,2X,I8)
40002 FORMAT(A,F12.8,2X,I8,2X,F12.8)
40003 FORMAT(5X,A,3(F20.15,2X))
40004 FORMAT(A,A,1X,A,1X,A)
50000 FORMAT(A,3(F18.12,2X),A)
50001 FORMAT(A,A,1X,A,1X,A,1X,A)
60000 FORMAT(A,A,2X,A,2X,A,3(2X,F10.5))

! end the timing
      CALL timestop(0.0E0_dp,handle)

!STOP

END SUBROUTINE mc_make_dat_file

END MODULE mc_misc

!******************************************************************************

