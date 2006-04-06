!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** mc_types.F *********************************************************
!!
!!   NAME
!!     mc_types
!!
!!   FUNCTION
!!     holds all the structure types needed for Monte Carlo, except
!!     the mc_environment_type
!!
!!   AUTHOR
!!     MJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE mc_types
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE fist_environment_types,          ONLY: fist_environment_type,&
                                             get_fist_env
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                             fist_nonbond_env_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: do_fist,&
                                             do_kg,&
                                             do_qs
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             molecule_kind_type
  USE pair_potential_types,            ONLY: pair_potential_pp_type
  USE physcon,                         ONLY: angstrom,&
                                             boltzmann,&
                                             joule
  USE string_utilities,                ONLY: xstring
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE 
! *** Global parameters ***

 PUBLIC :: mc_simulation_parameters_type,&
            mc_simulation_parameters_p_type,&
            mc_averages_type,mc_averages_p_type,&
            mc_moves_type, mc_moves_p_type,accattempt,&
            get_mc_par,set_mc_par,read_mc_section,&
            find_mc_rcut,determine_conf_prob


  TYPE mc_simulation_parameters_type
      PRIVATE
      INTEGER :: avbmc_atom
      INTEGER :: nstep
      INTEGER :: iupvolume
      INTEGER :: iuptrans
      INTEGER :: nbox
      INTEGER :: nchain
      INTEGER :: nunits
      INTEGER :: nunits_tot
      INTEGER :: nmoves
      INTEGER :: nvolmoves
      INTEGER :: nswapmoves
      INTEGER :: rm
      INTEGER :: cl
      INTEGER :: diff
      INTEGER :: nstart
      INTEGER :: source
      INTEGER :: group
      INTEGER :: iprint
      INTEGER :: print_level
      LOGICAL :: ldiscrete
      LOGICAL :: lbias
      LOGICAL :: ionode
      LOGICAL :: lrestart
      LOGICAL :: lstop
      CHARACTER ( LEN = 20 ) :: ensemble
      CHARACTER ( LEN = 40 ) :: restart_file_name
      CHARACTER ( LEN = 40 ) :: molecules_file
      CHARACTER ( LEN = 40 ) :: moves_file
      CHARACTER ( LEN = 40 ) :: coords_file
      CHARACTER ( LEN = 40 ) :: energy_file
      CHARACTER ( LEN = 40 ) :: displacement_file
      CHARACTER ( LEN = 40 ) :: cell_file
      CHARACTER ( LEN = 40 ) :: dat_file
      CHARACTER ( LEN = 40 ) :: data_file
      CHARACTER ( LEN = 200 ) :: box2_file
      CHARACTER ( LEN = 200 ) :: fft_lib
      CHARACTER ( LEN = 50 ) :: PROGRAM
      REAL ( dp ) :: discrete_step
      REAL ( dp ) :: rmvolume
      REAL ( dp ) :: rmbond
      REAL ( dp ) :: rmangle
      REAL ( dp ) :: rmdihedral
      REAL ( dp ) :: rmrot
      REAL ( dp ) :: rmtrans
      REAL ( dp ) :: temp
      REAL ( dp ) :: pressure
      REAL ( dp ) :: pbias
      REAL ( dp ) :: pmavbmc
      REAL ( dp ) :: pmswap
      REAL ( dp ) :: pmvolume
      REAL ( dp ) :: pmtraion
      REAL ( dp ) :: pmtrans
      REAL ( dp ) :: BETA
      REAL ( dp ) :: rcut
      REAL ( dp ) :: avbmc_rmin
      REAL ( dp ) :: avbmc_rmax
      REAL ( dp ),DIMENSION(1:3) :: conf_prob
      TYPE(section_vals_type),POINTER :: input_file
  END TYPE mc_simulation_parameters_type

  TYPE mc_simulation_parameters_p_type
      TYPE (mc_simulation_parameters_type),POINTER :: mc_par
  END TYPE mc_simulation_parameters_p_type

  TYPE mc_averages_type
      REAL(KIND = dp) :: ave_energy
      REAL(KIND = dp) :: ave_energy_squared
      REAL(KIND = dp) :: ave_volume
      REAL(KIND = dp) :: bond1
      REAL(KIND = dp) :: bond2
      REAL(KIND = dp) :: angle
      REAL(KIND = dp) :: molecules
  END TYPE mc_averages_type

  TYPE mc_averages_p_type
      TYPE (mc_averages_type),POINTER :: averages
  END TYPE mc_averages_p_type

  TYPE mc_moves_type
      TYPE ( accattempt ), POINTER :: bias_bond
      TYPE ( accattempt ), POINTER :: bias_angle
      TYPE ( accattempt ), POINTER :: bias_dihedral
      TYPE ( accattempt ), POINTER :: bias_trans
      TYPE ( accattempt ), POINTER :: bias_rot
      TYPE ( accattempt ), POINTER :: bond
      TYPE ( accattempt ), POINTER :: angle
      TYPE ( accattempt ), POINTER :: dihedral
      TYPE ( accattempt ), POINTER :: trans
      TYPE ( accattempt ), POINTER :: rot
      TYPE ( accattempt ), POINTER :: swap
      TYPE ( accattempt ), POINTER :: volume
      TYPE ( accattempt ), POINTER :: avbmc_inin
      TYPE ( accattempt ), POINTER :: avbmc_outin
      TYPE ( accattempt ), POINTER :: avbmc_inout
      TYPE ( accattempt ), POINTER :: avbmc_outout
      TYPE ( accattempt ), POINTER :: Quickstep
      REAL(KIND = dp) :: bias_energy,trans_dis,qtrans_dis
      INTEGER :: empty,grown,empty_conf,empty_avbmc
  END TYPE mc_moves_type

  TYPE accattempt
      INTEGER :: successes
      INTEGER :: qsuccesses
      INTEGER :: attempts
  END TYPE accattempt

  TYPE mc_moves_p_type
      TYPE(mc_moves_type), POINTER :: moves
  END TYPE mc_moves_p_type

! *** Global parameters ***
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mc_types'

CONTAINS
!****************************************************************************

 SUBROUTINE get_mc_par ( mc_par, nstep, iuptrans, iupvolume, nchain,&
       nunits,nunits_tot,nmoves,nvolmoves,nswapmoves,rm,cl,diff,nstart,&
       source,group,lbias,ionode,lrestart,lstop,rmvolume,rmbond,rmangle,&
       rmrot,rmtrans,temp,pressure,BETA,pmswap,pmvolume,pmtraion,pmtrans,&
       ensemble,PROGRAM,restart_file_name,molecules_file,moves_file,coords_file,&
       energy_file,displacement_file,cell_file,dat_file,data_file,box2_file,&
       fft_lib,iprint,rcut,ldiscrete,discrete_step,&
       conf_prob,pmavbmc,pbias,avbmc_atom,avbmc_rmin,avbmc_rmax,rmdihedral,&
       print_level,input_file)

! *****************************************************************************
!!****s* mc_types.F/get_mc_par
!!
!!   NAME
!!     get_mc_par
!!
!!   FUNCTION
!!     accesses the private elements of the mc_parameters_type
!!
!!   ARGUMENTS
!!     - mc_par: the structure mc parameters you want 
!!
!!   Suitable for parallel.
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    INTEGER, INTENT(OUT), OPTIONAL :: nstep, iuptrans, iupvolume, nchain, &
      nunits, nunits_tot, nmoves, nvolmoves, nswapmoves, rm, cl, diff, &
      nstart, source, group
    LOGICAL, INTENT(OUT), OPTIONAL           :: lbias, ionode, lrestart, lstop
    REAL(KIND=dp), INTENT(OUT), OPTIONAL :: rmvolume, rmbond, rmangle, rmrot, &
      rmtrans, temp, pressure, BETA, pmswap, pmvolume, pmtraion, pmtrans
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: ensemble, PROGRAM, &
      restart_file_name, molecules_file, moves_file, coords_file, &
      energy_file, displacement_file, cell_file, dat_file, data_file, &
      box2_file
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: fft_lib
    INTEGER, INTENT(OUT), OPTIONAL           :: iprint,print_level
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: rcut
    LOGICAL, INTENT(OUT), OPTIONAL           :: ldiscrete
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: discrete_step
    REAL(KIND=dp), DIMENSION(1:3), &
      INTENT(OUT), OPTIONAL                  :: conf_prob
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: pmavbmc, pbias
    INTEGER, INTENT(OUT), OPTIONAL           :: avbmc_atom
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: avbmc_rmin, avbmc_rmax, &
                                                rmdihedral
     TYPE(section_vals_type), POINTER, &
          OPTIONAL                           :: input_file

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_mc_par', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    IF ( PRESENT ( nstep ) ) nstep = mc_par % nstep
    IF ( PRESENT ( avbmc_atom ) ) avbmc_atom = mc_par % avbmc_atom
    IF ( PRESENT ( iuptrans ) ) iuptrans = mc_par % iuptrans
    IF ( PRESENT ( iupvolume ) ) iupvolume = mc_par % iupvolume
    IF ( PRESENT ( nchain ) ) nchain = mc_par % nchain
    IF ( PRESENT ( nunits ) ) nunits = mc_par % nunits
    IF ( PRESENT ( nunits_tot ) ) nunits_tot = mc_par % nunits_tot
    IF ( PRESENT ( nmoves ) ) nmoves = mc_par % nmoves
    IF ( PRESENT ( nvolmoves ) ) nvolmoves = mc_par % nvolmoves
    IF ( PRESENT ( nswapmoves ) ) nswapmoves = mc_par % nswapmoves
    IF ( PRESENT ( rm ) ) rm = mc_par % rm
    IF ( PRESENT ( cl ) ) cl = mc_par % cl
    IF ( PRESENT ( diff ) ) diff = mc_par % diff
    IF ( PRESENT ( nstart ) ) nstart = mc_par % nstart
    IF ( PRESENT ( source ) ) source = mc_par % source
    IF ( PRESENT ( group ) ) group = mc_par % group
    IF ( PRESENT ( iprint ) ) iprint = mc_par % iprint

    IF ( PRESENT ( lbias ) ) lbias = mc_par % lbias
    IF ( PRESENT ( ionode ) ) ionode = mc_par % ionode
    IF ( PRESENT ( lrestart ) ) lrestart = mc_par % lrestart
    IF ( PRESENT ( lstop ) ) lstop = mc_par % lstop
    IF ( PRESENT ( ldiscrete ) ) ldiscrete = mc_par % ldiscrete

    IF ( PRESENT ( conf_prob ) ) conf_prob(:) = mc_par % conf_prob(:)

    IF ( PRESENT ( avbmc_rmin ) ) avbmc_rmin = mc_par % avbmc_rmin
    IF ( PRESENT ( avbmc_rmax ) ) avbmc_rmax = mc_par % avbmc_rmax
    IF ( PRESENT ( rcut ) ) rcut = mc_par % rcut
    IF ( PRESENT ( discrete_step ) ) discrete_step = mc_par % discrete_step
    IF ( PRESENT ( rmvolume ) ) rmvolume = mc_par % rmvolume
    IF ( PRESENT ( rmbond ) ) rmbond = mc_par % rmbond
    IF ( PRESENT ( rmangle ) ) rmangle = mc_par % rmangle
    IF ( PRESENT ( rmdihedral ) ) rmdihedral = mc_par % rmdihedral
    IF ( PRESENT ( rmrot ) ) rmrot = mc_par % rmrot
    IF ( PRESENT ( rmtrans ) ) rmtrans = mc_par % rmtrans
    IF ( PRESENT ( temp ) ) temp = mc_par % temp
    IF ( PRESENT ( pressure ) ) pressure = mc_par % pressure
    IF ( PRESENT ( BETA ) ) BETA = mc_par % BETA
    IF ( PRESENT ( pmswap ) ) pmswap = mc_par % pmswap
    IF ( PRESENT ( pmvolume) ) pmvolume = mc_par % pmvolume
    IF ( PRESENT ( pmtraion ) ) pmtraion = mc_par % pmtraion
    IF ( PRESENT ( pmtrans ) ) pmtrans = mc_par % pmtrans
    IF ( PRESENT ( pmavbmc) ) pmavbmc = mc_par % pmavbmc
    IF ( PRESENT ( pbias) ) pbias = mc_par % pbias

    IF ( PRESENT ( ensemble ) ) ensemble = mc_par % ensemble
    IF ( PRESENT ( PROGRAM ) ) PROGRAM = mc_par % program
    IF ( PRESENT ( restart_file_name ) ) restart_file_name = &
      mc_par % restart_file_name
    IF ( PRESENT ( moves_file ) ) moves_file = mc_par % moves_file
    IF ( PRESENT ( coords_file ) ) coords_file = mc_par % coords_file
    IF ( PRESENT ( molecules_file ) ) molecules_file = mc_par % molecules_file
    IF ( PRESENT ( energy_file ) ) energy_file = mc_par % energy_file
    IF ( PRESENT ( displacement_file ) ) displacement_file = &
      mc_par % displacement_file
    IF ( PRESENT ( cell_file ) ) cell_file = mc_par % cell_file
    IF ( PRESENT ( dat_file ) ) dat_file = mc_par % dat_file
    IF ( PRESENT ( data_file ) ) data_file = mc_par % data_file
    IF ( PRESENT ( box2_file ) ) box2_file = mc_par % box2_file
    IF ( PRESENT ( fft_lib ) ) fft_lib = mc_par % fft_lib
    IF ( PRESENT ( print_level ) ) print_level = mc_par % print_level

    IF ( PRESENT ( input_file ) ) input_file => mc_par % input_file

  END SUBROUTINE get_mc_par

 SUBROUTINE set_mc_par ( mc_par, nchain,nunits,nunits_tot,rm,cl,&
      diff,nstart,rmvolume,rmbond,rmangle,rmdihedral,rmrot,rmtrans,PROGRAM,&
      nmoves,nswapmoves,lstop,temp,pressure,iuptrans,iupvolume,&
      pmswap,pmvolume,pmtraion,pmtrans,BETA,rcut,iprint,lbias,nstep,&
      lrestart,ldiscrete,discrete_step,pmavbmc)

! *****************************************************************************
!!****s* mc_types.F/set_mc_par
!!
!!   NAME
!!     set_mc_par
!!
!!   FUNCTION
!!     changes the private elements of the mc_parameters_type
!!
!!   ARGUMENTS
!!     - mc_par: the structure mc parameters you want 
!!
!!   Suitable for parallel.
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    INTEGER, INTENT(IN), OPTIONAL            :: nchain, nunits, nunits_tot, &
                                                rm, cl, diff, nstart
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: rmvolume, rmbond, rmangle, &
                                                rmdihedral, rmrot, rmtrans
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: PROGRAM
    INTEGER, INTENT(IN), OPTIONAL            :: nmoves, nswapmoves
    LOGICAL, INTENT(IN), OPTIONAL            :: lstop
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: temp, pressure
    INTEGER, INTENT(IN), OPTIONAL            :: iuptrans, iupvolume
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: pmswap, pmvolume, pmtraion, &
                                                pmtrans, BETA, rcut
    INTEGER, INTENT(IN), OPTIONAL            :: iprint
    LOGICAL, INTENT(IN), OPTIONAL            :: lbias
    INTEGER, INTENT(IN), OPTIONAL            :: nstep
    LOGICAL, INTENT(IN), OPTIONAL            :: lrestart, ldiscrete
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: discrete_step, pmavbmc

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_mc_par', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------
! These are the only values that change during the course of the simulation
! or are computed outside of this module

    IF ( PRESENT ( nstep ) ) mc_par % nstep = nstep
    IF ( PRESENT ( nchain ) ) mc_par % nchain = nchain
    IF ( PRESENT ( nunits ) ) mc_par % nunits = nunits
    IF ( PRESENT ( nunits_tot ) ) mc_par % nunits_tot = nunits_tot
    IF ( PRESENT ( rm ) ) mc_par % rm = rm
    IF ( PRESENT ( cl ) ) mc_par % cl = cl
    IF ( PRESENT ( diff ) ) mc_par % diff = diff
    IF ( PRESENT ( nstart ) ) mc_par % nstart = nstart
    IF ( PRESENT ( nmoves ) ) mc_par % nmoves = nmoves
    IF ( PRESENT ( nswapmoves ) ) mc_par % nswapmoves = nswapmoves
    IF ( PRESENT ( iprint ) ) mc_par % iprint = iprint
    IF ( PRESENT ( iuptrans ) ) mc_par % iuptrans = iuptrans
    IF ( PRESENT ( iupvolume ) ) mc_par % iupvolume = iupvolume

    IF ( PRESENT ( ldiscrete ) ) mc_par % ldiscrete = ldiscrete
    IF ( PRESENT ( lstop ) ) mc_par % lstop = lstop
    IF ( PRESENT ( lbias ) ) mc_par % lbias = lbias
    IF ( PRESENT ( lrestart ) ) mc_par % lrestart = lrestart

    IF ( PRESENT ( BETA ) ) mc_par % BETA = BETA
    IF ( PRESENT ( temp ) ) mc_par % temp = temp
    IF ( PRESENT ( rcut ) ) mc_par % rcut = rcut
    IF ( PRESENT ( pressure ) ) mc_par % pressure = pressure
    IF ( PRESENT ( pmvolume ) ) mc_par % pmvolume = pmvolume
    IF ( PRESENT ( pmswap ) ) mc_par % pmswap = pmswap
    IF ( PRESENT ( pmtrans ) ) mc_par % pmtrans = pmtrans
    IF ( PRESENT ( pmtraion ) ) mc_par % pmtraion = pmtraion
    IF ( PRESENT ( pmavbmc ) ) mc_par % pmavbmc = pmavbmc

    IF ( PRESENT ( discrete_step ) ) mc_par % discrete_step = discrete_step
    IF ( PRESENT ( rmvolume ) ) mc_par % rmvolume = rmvolume
    IF ( PRESENT ( rmbond ) ) mc_par % rmbond = rmbond
    IF ( PRESENT ( rmangle ) ) mc_par % rmangle = rmangle
    IF ( PRESENT ( rmdihedral ) ) mc_par % rmdihedral = rmdihedral
    IF ( PRESENT ( rmrot ) ) mc_par % rmrot = rmrot
    IF ( PRESENT ( rmtrans ) ) mc_par % rmtrans = rmtrans

    IF ( PRESENT ( PROGRAM ) ) mc_par % program = PROGRAM

  END SUBROUTINE set_mc_par

!******************************************************************************
!!****** mc/read_mc_section [1.0] *
!!
!!   NAME
!!     read_mc_section
!!
!!   FUNCTION
!!     reads in the Monte Carlo simulation parameters from an input file
!!
!!   ARGUMENTS
!!     - mc_par: the structure that will store the parameters
!!     - globenv: the global environment for the simulation
!!
!!   AUTHOR
!!     MJM
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &MC ... &END                                                     I
!!I                                                                           I
!!I program:     (traditional,gemc)       type of simulation
!!I nstep:                                number of steps
!!I lstop:       (yes,no)                 if no, peform nstep cycles
!!I ldiscrete:   (yes,no)                 if yes, perform a discrete volume
!!                                        move
!!I avbmc_rmin:                           the minimum distance used in AVBMC
!!I avbmc_rmax:                           the maximum distance used in AVBMC
!!I avbmc_atom:                           the atom used to calculate distances
!!                          between for the AVBMC swap move
!!I lbias:  (yes,no)                      pre-bias moves..requires an 
!!                                        additional file
!!I nmoves:                               number of moves between Quickstep
!!                                        calculations
!!I discrete_step:                        the size (in angstroms) of the cell
!!                                        length displacement for discrete
!!                                        volume_moves
!!I rmbond:                               maximum displacement for bond
!!                                        changes (in angstroms)
!!I rmangle:                              maximum displacement for angle
!!                                        changes (in degrees)
!!I rmdihedral:                           maximum displacement for dihedral
!!                                        changes (in degrees)
!!I rmtrans:                              maximum displacement for molecule
!!                                        translations (in angstroms)
!!I rmrot:                                maximum displacement for molecule
!!                                        rotations (in degrees)
!!I rmvolume:                             maximum displacement for volume
!!                                        changes (in cubic angstroms)
!!I iuptrans:                             every iuptrans step, update maximum
!!                                        displacements for angle,bond,rot,
!!                                        trans
!!I iupvolume:                            every iupvolume step, update maximum
!!                                        volume displacement
!!I conf_prob:                            holds the probabilities of doing a
!!                                        bond langth change (1), 
!!                                  angle change (2), or dihedral change(3)
!!I pmswap:                               probability of doing a swap move
!!I pmvolume:                             probability of doing a volume move
!!I pmavbmc:                              probability of doing an AVBMC1 move
!!I pmtraion:                             probability of doing a conformational
!!                                        move
!!I pbias:                                probability of doing different types
!!                                 of AVBMC moves (in-in,out-out,etc.)
!!                                        move
!!I pmtrans:                              probability of doing a translation
!!I temp:                                 temperature (in Kelvin)
!!I pressure:                             pressure (in bar)
!!I lrestart:     (yes,no)                read initial coordinates from a
!!                                        restart file
!!I restart_file_name:                    name of the restart file
!!I box2_file:                            name of the input file that contains
!!                                        the coords for the second box in
!!                                        a GEMC simulation
!!I---------------------------------------------------------------------------I
!!*****************************************************************************

SUBROUTINE read_mc_section ( mc_par, globenv, input_file_name, input_file,&
     error )


    TYPE(mc_simulation_parameters_type), &
      INTENT(OUT)                            :: mc_par
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(section_vals_type), POINTER         :: input_file
    CHARACTER(LEN=*),INTENT(IN)              :: input_file_name

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'read_mc_section'

    INTEGER                                  :: handle, ia, ie, iw, stop_num
    TYPE(section_vals_type), POINTER         :: mc_section
    LOGICAL                                  :: failure

!------------------------------------------------------------------------------
! begin the timing of the subroutine

    CPPrecondition(ASSOCIATED(input_file),cp_failure_level,routineN,error,failure)

    
      CALL timeset(routineN,'I','',handle)

      NULLIFY(mc_section)
      mc_section => section_vals_get_subs_vals(input_file,&
         "MOTION%MC",error=error)

! need the input file sturcutre that we're reading from for when we make
! dat files
      mc_par % input_file => input_file

! set the ionode and mepos
      mc_par % ionode = globenv % ionode
      mc_par % group = globenv % group
      mc_par % source = globenv % source

!..defaults
      mc_par % rmbond = 0.1d0
      mc_par % rmangle = 5.0d0
      mc_par % rmdihedral = 5.0d0
      mc_par % pbias = 0.5d0
      mc_par % pmavbmc = 0.0d0
      mc_par % pmswap = 0.0d0
      mc_par % pmvolume = 0.0d0
      mc_par % avbmc_atom = 1
      mc_par % avbmc_rmin = 2.0_dp
      mc_par % avbmc_rmax = 4.0_dp
      mc_par % nstart = 0
      mc_par % nvolmoves = 1
      mc_par % nswapmoves = 1
      mc_par % ldiscrete = .FALSE.
      mc_par % discrete_step = 0.05d0
      SELECT CASE (globenv%prog_name_id)
      CASE (do_fist)
         mc_par % iprint = 100 
      CASE (do_kg)
         mc_par % iprint = 10
      CASE (do_qs)
         mc_par % iprint = 1
      END SELECT
 
      iw = globenv % scr
      IF(globenv%ionode) WRITE ( iw, * )

!..filenames
      mc_par % program = input_file_name
      CALL xstring ( mc_par % program, ia, ie )
      mc_par%coords_file = mc_par % program(ia:ie) // '.coordinates'
      mc_par%molecules_file = mc_par % program(ia:ie) // '.molecules'
      mc_par%moves_file = mc_par % program(ia:ie) // '.moves'
      mc_par%energy_file = mc_par % program(ia:ie) // '.energy'
      mc_par%cell_file = mc_par % program(ia:ie) // '.cell'
      mc_par%displacement_file= mc_par % program(ia:ie)&
      // '.max_displacements'
      mc_par%data_file = mc_par % program(ia:ie) // '.data'
      stop_num=ie-3
      mc_par%dat_file = mc_par % program(ia:stop_num) // 'dat'
      WRITE(6,*) 'evnkx,v c ',mc_par % program

! set them into the input parameter structure as the new defaults
      CALL section_vals_val_set(mc_section,"COORDINATE_FILE_NAME",&
             c_val=mc_par%coords_file,error=error)
      CALL section_vals_val_set(mc_section,"DATA_FILE_NAME",&
             c_val=mc_par%data_file,error=error)
      CALL section_vals_val_set(mc_section,"CELL_FILE_NAME",&
             c_val=mc_par%cell_file,error=error)
      CALL section_vals_val_set(mc_section,"MAX_DISP_FILE_NAME",&
             c_val=mc_par%displacement_file,error=error)
      CALL section_vals_val_set(mc_section,"MOVES_FILE_NAME",&
             c_val=mc_par%moves_file,error=error)
      CALL section_vals_val_set(mc_section,"MOLECULES_FILE_NAME",&
             c_val=mc_par%molecules_file,error=error)
      CALL section_vals_val_set(mc_section,"ENERGY_FILE_NAME",&
             c_val=mc_par%energy_file,error=error)

! grab the FFT library name and print level...this is used for writing the dat file
! and hopefully will be changed
      mc_par % fft_lib = globenv % default_fft_library
      mc_par % print_level = globenv%print_level

 ! first, grab all the integer values
      CALL section_vals_val_get(mc_section,"NSTEP",i_val=mc_par%nstep,error=error)
      CALL section_vals_val_get(mc_section,"AVBMC_ATOM",i_val=mc_par%avbmc_atom,error=error)
      CALL section_vals_val_get(mc_section,"NMOVES",i_val=mc_par%nmoves,error=error)
      CALL section_vals_val_get(mc_section,"NVOLMOVES",i_val=mc_par%nvolmoves,error=error)
      CALL section_vals_val_get(mc_section,"NSWAPMOVES",i_val=mc_par%nswapmoves,error=error)
      CALL section_vals_val_get(mc_section,"IUPVOLUME",i_val=mc_par%iupvolume,error=error)
      CALL section_vals_val_get(mc_section,"IUPTRANS",i_val=mc_par%iuptrans,error=error)
      CALL section_vals_val_get(mc_section,"IPRINT",i_val=mc_par%iprint,error=error)

! now the real values
      CALL section_vals_val_get(mc_section,"AVBMC_RMIN",r_val=mc_par%avbmc_rmin,error=error)
      CALL section_vals_val_get(mc_section,"AVBMC_RMAX",r_val=mc_par%avbmc_rmax,error=error)
      CALL section_vals_val_get(mc_section,"PRESSURE",r_val=mc_par%pressure,error=error)
      CALL section_vals_val_get(mc_section,"TEMP",r_val=mc_par%temp,error=error)
      CALL section_vals_val_get(mc_section,"PMVOLUME",r_val=mc_par%pmvolume,error=error)
      CALL section_vals_val_get(mc_section,"PMSWAP",r_val=mc_par%pmswap,error=error)
      CALL section_vals_val_get(mc_section,"PMAVBMC",r_val=mc_par%pmavbmc,error=error)
      CALL section_vals_val_get(mc_section,"PBIAS",r_val=mc_par%pbias,error=error)
      CALL section_vals_val_get(mc_section,"PMTRANS",r_val=mc_par%pmtrans,error=error)
      CALL section_vals_val_get(mc_section,"PMTRAION",r_val=mc_par%pmtraion,error=error)
      CALL section_vals_val_get(mc_section,"DISCRETE_STEP",r_val=mc_par%discrete_step,error=error)
      CALL section_vals_val_get(mc_section,"RMVOLUME",r_val=mc_par%rmvolume,error=error)
      CALL section_vals_val_get(mc_section,"RMBOND",r_val=mc_par%rmbond,error=error)
      CALL section_vals_val_get(mc_section,"RMANGLE",r_val=mc_par%rmangle,error=error)
      CALL section_vals_val_get(mc_section,"RMDIHEDRAL",r_val=mc_par%rmdihedral,error=error)
      CALL section_vals_val_get(mc_section,"RMROT",r_val=mc_par%rmrot,error=error)
      CALL section_vals_val_get(mc_section,"RMTRANS",r_val=mc_par%rmtrans,error=error)


! now some logical values
      CALL section_vals_val_get(mc_section,"LBIAS",l_val=mc_par%lbias,error=error)
      CALL section_vals_val_get(mc_section,"LDISCRETE",l_val=mc_par%ldiscrete,error=error)
      CALL section_vals_val_get(mc_section,"LSTOP",l_val=mc_par%lstop,error=error)
      CALL section_vals_val_get(mc_section,"RESTART",l_val=mc_par%lrestart,error=error)
      CALL section_vals_val_get(mc_section,"LBIAS",l_val=mc_par%lbias,error=error)

! finally the character values
      CALL section_vals_val_get(mc_section,"PROGRAM",c_val=mc_par%ensemble,error=error)
      CALL section_vals_val_get(mc_section,"RESTART_FILE_NAME",c_val=mc_par%restart_file_name,error=error)
      CALL section_vals_val_get(mc_section,"COORDINATE_FILE_NAME",c_val=mc_par%coords_file,error=error)
      CALL section_vals_val_get(mc_section,"ENERGY_FILE_NAME",c_val=mc_par%energy_file,error=error)
      CALL section_vals_val_get(mc_section,"MOVES_FILE_NAME",c_val=mc_par%moves_file,error=error)
      CALL section_vals_val_get(mc_section,"MOLECULES_FILE_NAME",c_val=mc_par%molecules_file,error=error)
      CALL section_vals_val_get(mc_section,"CELL_FILE_NAME",c_val=mc_par%cell_file,error=error)
      CALL section_vals_val_get(mc_section,"DATA_FILE_NAME",c_val=mc_par%data_file,error=error)
      CALL section_vals_val_get(mc_section,"MAX_DISP_FILE_NAME",c_val=mc_par%displacement_file,error=error)
      CALL section_vals_val_get(mc_section,"BOX2_FILE_NAME",c_val=mc_par%box2_file,error=error)

!..end of parsing the input section

!..write some information to output
      IF (mc_par % ionode ) THEN
         WRITE ( iw, '( A )' ) ' MC| Monte Carlo Protocol '
         WRITE ( iw, '( A,T71,I10 )' ) ' MC| total number of steps ', &
            mc_par % nstep
         WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| pmvolume ', &
            mc_par % pmvolume
         WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| pmswap ', &
            mc_par % pmswap
         WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| pmavbmc ', &
            mc_par % pmavbmc
         WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| pbias ', &
            mc_par % pbias
         WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| pmtraion ', &
            mc_par % pmtraion
         WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| pmtrans ', &
            mc_par % pmtrans
         WRITE ( iw, '( A,T71,I10 )' ) ' MC| iupvolume ', &
            mc_par % iupvolume
         WRITE ( iw, '( A,T71,I10 )' ) ' MC| iuptrans ', &
            mc_par % iuptrans
         WRITE ( iw, '( A,T71,I10 )' ) ' MC| iprint ', &
            mc_par % iprint
         WRITE ( iw, '( A,T58,A20 )' ) ' MC| ensemble ', &
            mc_par % ensemble
         IF (mc_par%ensemble .EQ. 'GEMC-NVT') THEN
            WRITE ( iw, '( A,T58,A)' ) ' MC| Box 2 file', &
            TRIM(mc_par % box2_file)
         ENDIF
         WRITE ( iw, '( A,T58,A )' ) ' MC| Name of restart file:',&
            TRIM(mc_par % restart_file_name)
         WRITE ( iw, '( A,A,T44,A )' ) ' MC| Name of output ',&
            'coordinate file:',&
         TRIM(mc_par % coords_file)
         WRITE ( iw, '( A,T44,A )' ) ' MC| Name of output data file:',&
            TRIM(mc_par % data_file)
         WRITE ( iw, '( A,A,T44,A )' ) ' MC| Name of output ',&
            'molecules file:',&
            TRIM(mc_par %molecules_file)
         WRITE ( iw, '( A,T44,A )' ) ' MC| Name of output moves file:',&
            TRIM(mc_par % moves_file)
         WRITE ( iw, '( A,T44,A )') ' MC| Name of output energy file:',&
            TRIM(mc_par % energy_file)
         WRITE ( iw, '( A,T44,A )' ) ' MC| Name of output cell file:',&
            TRIM(mc_par % cell_file)
         WRITE ( iw, '( A,A,T44,A )' ) ' MC| Name of output',&
            ' displacement file:',&
            TRIM(mc_par % displacement_file)
         IF(mc_par % ldiscrete) THEN
            WRITE ( iw, '(A,A,T71,F10.3)' ) ' MC| discrete step size',&
            '[angstroms]', &
            mc_par % discrete_step
         ELSE
            WRITE ( iw, '( A,A,T71,F10.3 )' ) ' MC| rmvolume ',&
            '[cubic angstroms]', &
            mc_par % rmvolume
         ENDIF
         WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| rmbond [angstroms]', &
            mc_par % rmbond
         WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| rmangle [degrees]', &
            mc_par % rmangle
         WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| rmdihedral [degrees]', &
            mc_par % rmdihedral
         WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| rmtrans [angstroms]', &
            mc_par % rmtrans
         WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| rmrot [degrees]', &
            mc_par % rmrot
         WRITE ( iw, '( A,T71,F10.2 )' ) ' MC| Temperature [K] ', &
            mc_par % temp
         WRITE ( iw, '( A,T71,F10.5 )' ) ' MC| Pressure [bar] ', &
            mc_par % pressure
         WRITE ( iw, '( A,T71,I10 )' ) ' MC| AVBMC target atom ', &
            mc_par % avbmc_atom
         WRITE ( iw, '( A,T71,F10.5 )' ) ' MC| AVBMC inner cutoff [ang]', &
            mc_par % avbmc_rmin
         WRITE ( iw, '( A,T71,F10.5 )' ) ' MC| AVBMC outer cutoff [ang]', &
            mc_par % avbmc_rmax
         IF ( mc_par % lrestart ) THEN
            WRITE ( iw, '(A,A)') ' MC| Initial data will be read from a',&
            ' restart file.'
         ENDIF
         IF ( mc_par % lbias ) THEN
            WRITE ( iw, '(A)') ' MC| The moves will be biased.'
         ELSE
            WRITE ( iw, '(A)') ' MC| The moves will not be biased,'
         ENDIF
         IF (mc_par%nmoves .EQ. 1) THEN
            WRITE(iw,'(A,A)') ' MC| A full energy calculation ',&
            'will be done at every step.'
         ELSE
            WRITE( iw, '(A,I4,A,A)' ) ' MC| ',mc_par%nmoves,&
            ' moves will be attempted ',&
            'before a Quickstep energy calculation'
            WRITE( iw, '(A)' ) ' MC|      takes place.'
         ENDIF
         IF (mc_par%pmvolume .GT. 0.0e0_dp ) THEN
            WRITE( iw, '(A,I4,A,A)') ' MC| ',mc_par%nvolmoves,&
            ' volume moves will be attempted ',&
            'before a Quickstep energy calculation'
            WRITE( iw, '(A)' ) ' MC|      takes place.'
         ENDIF
         IF (mc_par%pmswap .GT. 0.0e0_dp ) THEN
            WRITE( iw, '(A,I4,A,A)') ' MC| ',mc_par%nswapmoves,&
            ' classical swap insertions will be attempted ',&
            'per Quickstep swap move'
         ENDIF
      END IF

! figure out what beta (1/kT) is in atomic units (1/Hartree)
      mc_par % BETA = 1 / mc_par%temp / boltzmann * joule

! convert from degrees to radians
      mc_par%rmrot = mc_par%rmrot/180.0e0_dp*pi
      mc_par%rmangle = mc_par%rmangle/180.0e0_dp*pi
      mc_par%rmdihedral = mc_par%rmdihedral/180.0e0_dp*pi
! convert from bar to a.u.
      CALL convert_to_cp2k_units(pressure=mc_par%pressure)
! convert from angstrom to a.u.
      CALL convert_to_cp2k_units(unit_char='ANGSTROM',&
         length=mc_par%rmtrans)
      CALL convert_to_cp2k_units(unit_char='ANGSTROM',&
         length=mc_par%rmbond)
      CALL convert_to_cp2k_units(unit_char='ANGSTROM',&
         length=mc_par%rmvolume,l_power=3)
      CALL convert_to_cp2k_units(unit_char='ANGSTROM',&
         length=mc_par%discrete_step)
      CALL convert_to_cp2k_units(unit_char='ANGSTROM',&
         length=mc_par%avbmc_rmin)
      CALL convert_to_cp2k_units(unit_char='ANGSTROM',&
         length=mc_par%avbmc_rmax)

! end the timing
      CALL timestop(0.0e0_dp,handle)


END SUBROUTINE read_mc_section

!******************************************************************************
!!****** mc/find_mc_rcut [1.0] *
!!
!!   NAME
!!     find_mc_rcut
!!
!!   FUNCTION
!!     finds the largest interaction cutoff value in a classical simulation
!!     so we know the smallest size we can make the box in a volume move
!!
!!   ARGUMENTS
!!     - mc_par: the structure that will store the parameters
!!     - force_env: the force environment that we'll grab the rcut parameter
!!                  out of
!!     - lterminate: set to .TRUE. if one of the sides of the box is
!!           less than twice the cutoff
!!
!!   Suitable for parallel.
!!
!!   AUTHOR
!!     MJM
!!
!!   NOTES
!!*****************************************************************************

SUBROUTINE find_mc_rcut ( mc_par, force_env, lterminate )


    TYPE(mc_simulation_parameters_type), &
      INTENT(INOUT)                          :: mc_par
    TYPE(force_env_type), POINTER            :: force_env
    LOGICAL, INTENT(OUT)                     :: lterminate

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'find_mc_rcut'

    INTEGER                                  :: itype, jtype
    REAL(KIND=dp)                            :: rcutsq_max
    REAL(KIND=dp), DIMENSION(1:3)            :: abc
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type)                      :: error
    TYPE(fist_environment_type), POINTER     :: fist_env
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(pair_potential_pp_type), POINTER    :: potparm

!------------------------------------------------------------------------------

      NULLIFY(cell,potparm,fist_nonbond_env,fist_env)

      lterminate=.FALSE.
      CALL force_env_get(force_env,fist_env=fist_env)
      CALL get_fist_env(fist_env,cell=cell,&
         fist_nonbond_env=fist_nonbond_env,error=error)
      CALL fist_nonbond_env_get (fist_nonbond_env, potparm=potparm)
      CALL get_cell(cell,abc=abc)
      
! find the largest value of rcutsq
      rcutsq_max=0.0d0
      DO itype=1,SIZE(potparm%pot,1)
         DO jtype=itype,SIZE(potparm%pot,2)
            IF(potparm%pot(itype,jtype)%pot%rcutsq .GT. rcutsq_max) &
                rcutsq_max=potparm%pot(itype,jtype)%pot%rcutsq
         ENDDO
      ENDDO

! check to make sure all box dimensions are greater than two times this
! value
      mc_par % rcut=rcutsq_max**0.5_dp
      DO itype=1,3
         IF(abc(itype) .LT. 2.0_dp*mc_par % rcut) THEN
            lterminate=.TRUE.
         ENDIF
      ENDDO

END SUBROUTINE find_mc_rcut

!******************************************************************************
!!****** mc/determine_conf_prob [1.0] *
!!
!!   NAME
!!     determine_conf_prob
!!
!!   FUNCTION
!!     determines the probabilities of doing interal conformational moves 
!!     in the system, based on how many of each type there are...the 
!!     result is stored in mc_par % conf_prob
!!
!!   ARGUMENTS
!!     - mc_par: the structure that will store the parameters
!!     - molecule_kind: the pointer to the molecule kind we're looking at
!!
!!   Suitable for parallel.
!!
!!   AUTHOR
!!     MJM
!!
!!   NOTES
!!*****************************************************************************

SUBROUTINE determine_conf_prob ( mc_par, molecule_kind)


    TYPE(mc_simulation_parameters_type), &
      INTENT(INOUT)                          :: mc_par
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    CHARACTER(LEN=*), PARAMETER :: routineN = 'determine_conf_prob'

    INTEGER                                  :: nbend, nbond, ntorsion, total

!------------------------------------------------------------------------------

      CALL get_molecule_kind(molecule_kind,nbond=nbond,nbend=nbend,&
         ntorsion=ntorsion)

      total=nbond+nbend+ntorsion

      mc_par % conf_prob(1)=REAL(nbond,dp)/REAL(total,dp)
      mc_par % conf_prob(2)=REAL(nbend,dp)/REAL(total,dp)
      mc_par % conf_prob(3)=REAL(ntorsion,dp)/REAL(total,dp)

END SUBROUTINE determine_conf_prob

END MODULE mc_types

!******************************************************************************
