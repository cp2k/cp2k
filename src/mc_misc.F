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
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
!  USE cp_control_types,                ONLY: dft_control_type
!  USE cp_control_utils,                ONLY: xc_get_rho_smooth_label,&
!                                             xc_get_deriv_method_label
  USE external_potential_types,        ONLY: get_potential,&
                                             gth_potential_type,&
                                             elp_potential_type
  USE ewald_environment_types ,        ONLY: ewald_environment_type,&
                                             ewald_env_get
  USE fist_environment_types,          ONLY: fist_environment_type,&
                                             get_fist_env
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_type,&
                                             fist_nonbond_env_get
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type,&
                                             use_fist_force,&
                                             use_kg_force,&
                                             use_qs_force
  USE input_constants,                 ONLY: do_ewald_none,&
                                             do_ewald_ewald,&
                                             do_ewald_pme,&
                                             do_ewald_spme
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_get_subs_vals2,&
                                             section_vals_get_subs_vals
  USE kg_environment_types,            ONLY: kg_environment_type,&
                                             get_kg_env
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: pi
  USE mc_types,                        ONLY: mc_averages_type,&
                                             mc_moves_type,&
                                             mc_simulation_parameters_type,&
                                             get_mc_par
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             molecule_kind_type,&
                                             atom_type
  USE pair_potential_types,            ONLY: pair_potential_pp_type
  USE particle_list_types,             ONLY: particle_list_type
  USE physcon,                         ONLY: angstrom,&
                                             n_avogadro,&
                                             kelvin
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cell_types,                 ONLY: cell_type,&
                                             get_cell
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

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
!!   Suitable for parallel.
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

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "final_mc_write"

    LOGICAL                                  :: lbias
    INTEGER                                  :: handle,nvolmoves
    REAL(KIND=dp)                            :: pmswap,rmvolume,rmbond,&
      rmangle,rmdihedral,rmtrans,rmrot
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
      force_env,fft_lib)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: coordinates
    INTEGER, INTENT(IN)                      :: natoms_tot
    REAL(KIND=dp), DIMENSION(1:3),INTENT(IN) :: box_length
    CHARACTER(LEN=*), INTENT(IN)             :: filename
    CHARACTER(LEN=*), INTENT(IN)             :: fft_lib
    TYPE(force_env_type),POINTER             :: force_env

    TYPE(cp_error_type)                      :: error
    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_make_dat_file", &
      routineP = module_name//':'//routineN

    TYPE(cell_type),POINTER                  :: cell_ref,cell
    TYPE(qs_environment_type),POINTER        :: qs_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(gth_potential_type), POINTER      :: gth_potential
    TYPE(elp_potential_type), POINTER      :: elp_potential
    TYPE(pair_potential_pp_type), POINTER      :: potparm
    TYPE(fist_nonbond_env_type), POINTER      :: fist_nonbond_env
    TYPE(gto_basis_set_type), POINTER      :: orb_basis_set
    TYPE ( kg_environment_type ), &
      POINTER                                :: kg_env                     
    TYPE ( fist_environment_type ), &
      POINTER                                :: fist_env 
    TYPE ( ewald_environment_type ), POINTER      :: ewald_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(atom_type), DIMENSION(:),POINTER    :: atom_list
  
    INTEGER                                  :: handle, iatom, unit,&
      itype,ntypes,gmax,o_spline,ns_max,ewald_type,jtype,ichain,&
      nchains,natoms,i_val,ifunct
    LOGICAL                                  :: use_ref_cell,failure
    REAL(KIND=dp),DIMENSION(1:3)             :: abc
    REAL(KIND=dp)                            :: alpha,epsilon,charge,r_val
    CHARACTER(LEN=5)                         :: cfile
    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=default_string_length)     :: name,unit_type,c_val
    TYPE(section_vals_type), POINTER         :: input,xc_fun,xc_fun_section

!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! nullify some stuff
      NULLIFY(qs_env,atomic_kind_set,atomic_kind,gth_potential,&
         orb_basis_set,elp_potential,input,xc_fun,kg_env)

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
         natom=natoms)
      
! find out how many atom types we have
      ntypes=SIZE(atomic_kind_set(:))

! open the file
      CALL open_file(file_name=filename,unit_number=unit,&
         file_action='WRITE',file_status='REPLACE')

! write all the stuff down
      WRITE(unit,10000) "&FORCE_EVAL"
      SELECT CASE ( force_env%in_use )
      CASE ( use_fist_force )
! grab some of the Ewald keywords
         WRITE(unit,10000) "  &MM"
         WRITE(unit,10000) "    &POISSON_MM"
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
         WRITE(unit,10000) "    &END POISSON_MM"


! rcut
         CALL fist_nonbond_env_get (fist_nonbond_env, potparm=potparm)

! now some of the force field stuff
         WRITE(unit,10000) "    &FORCEFIELD"
! doesn't look like these values are stored past read-in
         WRITE(unit,10000) "      parmfile water.pot"
         WRITE(unit,10000) "      parmtype CHM"

         
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
         DO ifunct=1,2
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
!         WRITE(unit,10000) "      EXTRAPOLATION USE_PREV_WF"
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
         WRITE(unit,10000) "      &XC_FUNCTIONAL"
         DO ifunct=1,2
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
      WRITE(unit,10000) '      CONN_FILE topology_' &
      // TRIM(ADJUSTL(cfile)) // '.psf'
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
      WRITE(unit,20001) ' FFTLIB       ',fft_lib
      WRITE(unit,10000) ' PRINT_LEVEL  LOW'
      WRITE(unit,10000) "&END"

      SELECT CASE ( force_env%in_use )
      CASE ( use_fist_force )

      CASE ( use_qs_force )

! let's try to deal with the reference cell
         IF(use_ref_cell) THEN
            WRITE(unit,10000) '&CELL_REF'
            WRITE(unit,10000) ' UNIT         ANGSTROM'
            WRITE(unit,40000) '  ABC          ',&
               abc(1:3)*angstrom
            WRITE(unit,10000) "&END"
            WRITE(unit,*) 
         ENDIF
      CASE ( use_kg_force )
! let's try to deal with the reference cell
         IF(use_ref_cell) THEN
            WRITE(unit,10000) '&CELL_REF'
            WRITE(unit,10000) ' UNIT         ANGSTROM'
            WRITE(unit,40000) '  ABC          ',&
               abc(1:3)*angstrom
            WRITE(unit,10000) "&END"
            WRITE(unit,*) 
         ENDIF
      CASE default
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
      END SELECT

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
50000 FORMAT(A,3(F18.12,2X),A)
60000 FORMAT(A,A,2X,A,2X,A,3(2X,F10.5))

! end the timing
      CALL timestop(0.0E0_dp,handle)

!STOP

END SUBROUTINE mc_make_dat_file

END MODULE mc_misc

!******************************************************************************

