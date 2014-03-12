!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief holds all the structure types needed for Monte Carlo, except
!>      the mc_environment_type
!> \par History
!>      none
!> \author MJM
! *****************************************************************************
MODULE mc_types
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE f77_blas
  USE fist_environment_types,          ONLY: fist_env_get,&
                                             fist_environment_type
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                             fist_nonbond_env_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_type,&
                                             use_fist_force,&
                                             use_qs_force
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: do_fist,&
                                             do_qs
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp,&
                                             dp_size,&
                                             int_size
  USE mathconstants,                   ONLY: pi
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_p_type
  USE molecule_kind_types,             ONLY: atom_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE pair_potential_types,            ONLY: pair_potential_pp_type
  USE particle_list_types,             ONLY: particle_list_p_type
  USE physcon,                         ONLY: boltzmann,&
                                             joule
  USE string_utilities,                ONLY: uppercase,&
                                             xstring
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

 PUBLIC :: mc_simpar_type,&
            mc_simulation_parameters_p_type,&
            mc_averages_type,mc_averages_p_type,&
            mc_moves_type, mc_moves_p_type,accattempt,&
            get_mc_par,set_mc_par,read_mc_section,&
            find_mc_rcut,mc_molecule_info_type,&
            mc_determine_molecule_info,mc_molecule_info_destroy,&
            mc_sim_par_create,mc_sim_par_destroy,&
            get_mc_molecule_info,&
            mc_input_file_type,mc_input_file_create,&
            get_mc_input_file,mc_input_file_destroy,&
            mc_input_parameters_check,&
            mc_ekin_type

! *****************************************************************************
  TYPE mc_simpar_type
! contains all the parameters needed for running Monte Carlo simulations
      PRIVATE
      INTEGER,DIMENSION(:),POINTER :: avbmc_atom
      INTEGER :: nstep
      INTEGER :: iupvolume
      INTEGER :: iuptrans
      INTEGER :: nvirial
      INTEGER :: nbox
      INTEGER :: nmoves
      INTEGER :: nswapmoves
      INTEGER :: rm
      INTEGER :: cl
      INTEGER :: diff
      INTEGER :: nstart
      INTEGER :: source
      INTEGER :: group
      INTEGER :: iprint
      LOGICAL :: ldiscrete
      LOGICAL :: lbias
      LOGICAL :: ionode
      LOGICAL :: lrestart
      LOGICAL :: lstop
      LOGICAL :: lhmc
      CHARACTER ( LEN = 20 ) :: ensemble
      CHARACTER ( LEN = default_path_length ) :: restart_file_name
      CHARACTER ( LEN = default_path_length ) :: molecules_file
      CHARACTER ( LEN = default_path_length ) :: moves_file
      CHARACTER ( LEN = default_path_length ) :: coords_file
      CHARACTER ( LEN = default_path_length ) :: energy_file
      CHARACTER ( LEN = default_path_length ) :: displacement_file
      CHARACTER ( LEN = default_path_length ) :: cell_file
      CHARACTER ( LEN = default_path_length ) :: dat_file
      CHARACTER ( LEN = default_path_length ) :: data_file
      CHARACTER ( LEN = default_path_length ) :: box2_file
      CHARACTER ( LEN = 200 ) :: fft_lib
      CHARACTER ( LEN = 50 ) :: PROGRAM
      REAL ( dp ),DIMENSION(:),POINTER :: pmtraion_mol,pmtrans_mol,pmrot_mol,pmswap_mol,&
           pbias,pmavbmc_mol
      REAL ( dp ) :: discrete_step
      REAL ( dp ) :: rmvolume
      REAL ( dp ),DIMENSION(:),POINTER :: rmbond,rmangle,rmdihedral,rmrot,rmtrans
      REAL ( dp ),DIMENSION(:),POINTER :: eta
      REAL ( dp ) :: temperature
      REAL ( dp ) :: pressure
      REAL ( dp ) :: pmavbmc
      REAL ( dp ) :: pmswap
      REAL ( dp ) :: pmvolume,pmvol_box
      REAL ( dp ) :: pmhmc,pmhmc_box
      REAL ( dp ) :: pmtraion
      REAL ( dp ) :: pmtrans
      REAL ( dp ) :: BETA
      REAL ( dp ) :: rcut
      REAL ( dp ),DIMENSION(:),POINTER :: avbmc_rmin,avbmc_rmax
      REAL ( dp ),DIMENSION(:),POINTER :: virial_temps
      TYPE(mc_input_file_type),POINTER ::&
           mc_input_file,mc_bias_file
      TYPE(section_vals_type),POINTER :: input_file
      TYPE(mc_molecule_info_type),POINTER :: mc_molecule_info
      REAL ( dp ) :: exp_min_val, exp_max_val, min_val, max_val
      INTEGER     :: rand2skip
  END TYPE mc_simpar_type

! *****************************************************************************
  TYPE mc_ekin_type
! contains the kinetic energy of an MD sequence from hybrid MC
     REAL ( dp ) :: initial_ekin,final_ekin
   END TYPE mc_ekin_type
! *****************************************************************************
  TYPE mc_input_file_type
! contains all the text of the input file
      PRIVATE
      INTEGER :: run_type_row,run_type_column,coord_row_start,coord_row_end,&
           cell_row,cell_column,global_row_end,in_use,nunits_empty,&
           motion_row_end,motion_row_start
      INTEGER,DIMENSION(:),POINTER :: mol_set_nmol_row,mol_set_nmol_column
      CHARACTER ( LEN = default_string_length ),DIMENSION(:),POINTER :: text
      CHARACTER ( LEN = default_string_length ),DIMENSION(:),POINTER :: atom_names_empty
      REAL(dp),DIMENSION(:,:),POINTER :: coordinates_empty
   END TYPE mc_input_file_type

! *****************************************************************************
  TYPE mc_molecule_info_type
! contains information on molecules...I had to do this because I use multiple force
! environments, and they know nothing about each other
      PRIVATE
      INTEGER :: nboxes
      CHARACTER ( LEN = default_string_length ),DIMENSION(:),POINTER :: names
      CHARACTER ( LEN = default_string_length ),&
           DIMENSION(:,:),POINTER :: atom_names
      REAL ( dp ),DIMENSION(:,:),POINTER :: conf_prob,mass
      INTEGER, DIMENSION(:,:),POINTER :: nchains
      INTEGER :: nmol_types,nchain_total
      INTEGER,DIMENSION(:),POINTER :: nunits,mol_type,nunits_tot,in_box
 END TYPE mc_molecule_info_type

! *****************************************************************************
  TYPE mc_simulation_parameters_p_type
! a pointer for multiple boxes
      TYPE (mc_simpar_type),POINTER :: mc_par
  END TYPE mc_simulation_parameters_p_type

! *****************************************************************************
  TYPE mc_averages_type
! contains some data that is averaged throughout the course of a run
      REAL(KIND = dp) :: ave_energy
      REAL(KIND = dp) :: ave_energy_squared
      REAL(KIND = dp) :: ave_volume
      REAL(KIND = dp) :: molecules
  END TYPE mc_averages_type

! *****************************************************************************
  TYPE mc_averages_p_type
      TYPE (mc_averages_type),POINTER :: averages
  END TYPE mc_averages_p_type

! *****************************************************************************
  TYPE mc_moves_type
! contains data for how many moves of a give type we've accepted/rejected
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
      TYPE ( accattempt ), POINTER :: hmc
      TYPE ( accattempt ), POINTER :: avbmc_inin
      TYPE ( accattempt ), POINTER :: avbmc_outin
      TYPE ( accattempt ), POINTER :: avbmc_inout
      TYPE ( accattempt ), POINTER :: avbmc_outout
      TYPE ( accattempt ), POINTER :: Quickstep
      REAL(KIND = dp) :: trans_dis,qtrans_dis
      INTEGER :: empty,grown,empty_conf,empty_avbmc
  END TYPE mc_moves_type

! *****************************************************************************
  TYPE accattempt
      INTEGER :: successes
      INTEGER :: qsuccesses
      INTEGER :: attempts
  END TYPE accattempt

! *****************************************************************************
  TYPE mc_moves_p_type
      TYPE(mc_moves_type), POINTER :: moves
  END TYPE mc_moves_p_type

! *** Global parameters ***
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mc_types'

CONTAINS

! *****************************************************************************
!> \brief accesses the private elements of the mc_input_file_type
!> \param mc_input_file the input file you want data for
!>
!>    Suitable for parallel.
!> \author MJM
! *****************************************************************************
 SUBROUTINE get_mc_input_file ( mc_input_file,run_type_row,run_type_column,&
      coord_row_start,coord_row_end,cell_row,cell_column,global_row_end,&
      mol_set_nmol_row,mol_set_nmol_column,in_use,text,atom_names_empty,&
      nunits_empty,coordinates_empty,motion_row_end,motion_row_start)

    TYPE(mc_input_file_type), POINTER        :: mc_input_file
    INTEGER, INTENT(OUT), OPTIONAL :: run_type_row, run_type_column, &
      coord_row_start, coord_row_end, cell_row, cell_column, global_row_end
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: mol_set_nmol_row, &
                                                mol_set_nmol_column
    INTEGER, INTENT(OUT), OPTIONAL           :: in_use
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), OPTIONAL, POINTER        :: text, atom_names_empty
    INTEGER, INTENT(OUT), OPTIONAL           :: nunits_empty
    REAL(dp), DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: coordinates_empty
    INTEGER, INTENT(OUT), OPTIONAL           :: motion_row_end, &
                                                motion_row_start

   IF(PRESENT(in_use)) in_use=mc_input_file%in_use
   IF(PRESENT(run_type_row)) run_type_row=mc_input_file%run_type_row
   IF(PRESENT(run_type_column)) run_type_column=mc_input_file%run_type_column
   IF(PRESENT(coord_row_start)) coord_row_start=mc_input_file%coord_row_start
   IF(PRESENT(coord_row_end)) coord_row_end=mc_input_file%coord_row_end
   IF(PRESENT(cell_row)) cell_row=mc_input_file%cell_row
   IF(PRESENT(cell_column)) cell_column=mc_input_file%cell_column
   IF(PRESENT(global_row_end)) global_row_end=mc_input_file%global_row_end
   IF(PRESENT(nunits_empty)) nunits_empty=mc_input_file%nunits_empty
   IF(PRESENT(motion_row_start)) motion_row_start=mc_input_file%motion_row_start
   IF(PRESENT(motion_row_end)) motion_row_end=mc_input_file%motion_row_end

   IF(PRESENT(mol_set_nmol_row)) mol_set_nmol_row=> mc_input_file%mol_set_nmol_row
   IF(PRESENT(mol_set_nmol_column)) mol_set_nmol_column=> mc_input_file%mol_set_nmol_column
   IF(PRESENT(text)) text=> mc_input_file%text
   IF(PRESENT(atom_names_empty)) atom_names_empty=> mc_input_file%atom_names_empty
   IF(PRESENT(coordinates_empty)) coordinates_empty=> mc_input_file%coordinates_empty

 END SUBROUTINE GET_MC_INPUT_FILE
! *****************************************************************************
 SUBROUTINE get_mc_par ( mc_par, nstep, nvirial, iuptrans, iupvolume, &
       nmoves,nswapmoves,rm,cl,diff,nstart,&
       source,group,lbias,ionode,lrestart,lstop,rmvolume,rmbond,rmangle,&
       rmrot,rmtrans,temperature,pressure,BETA,pmswap,pmvolume,pmtraion,pmtrans,&
       ensemble,PROGRAM,restart_file_name,molecules_file,moves_file,coords_file,&
       energy_file,displacement_file,cell_file,dat_file,data_file,box2_file,&
       fft_lib,iprint,rcut,ldiscrete,discrete_step,&
       pmavbmc,pbias,avbmc_atom,avbmc_rmin,avbmc_rmax,rmdihedral,&
       input_file,mc_molecule_info,pmswap_mol,pmavbmc_mol,pmtrans_mol,pmrot_mol,&
       pmtraion_mol,mc_input_file,mc_bias_file,pmvol_box,virial_temps,&
       exp_min_val, exp_max_val, min_val, max_val, eta, pmhmc, pmhmc_box,lhmc,rand2skip)


! *****************************************************************************
!> \brief accesses the private elements of the mc_parameters_type
!> \param mc_par the structure mc parameters you want
!>
!>    Suitable for parallel.
!> \author MJM
    TYPE(mc_simpar_type), POINTER            :: mc_par
    INTEGER, INTENT(OUT), OPTIONAL           :: nstep, nvirial, iuptrans, &
                                                iupvolume, nmoves, &
                                                nswapmoves, rm, cl, diff, &
                                                nstart, source, group
    LOGICAL, INTENT(OUT), OPTIONAL           :: lbias, ionode, lrestart, lstop
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: rmvolume
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: rmbond, rmangle, rmrot, &
                                                rmtrans
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: temperature, pressure, BETA, &
                                                pmswap, pmvolume, pmtraion, &
                                                pmtrans
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: ensemble, PROGRAM, &
      restart_file_name, molecules_file, moves_file, coords_file, &
      energy_file, displacement_file, cell_file, dat_file, data_file, &
      box2_file, fft_lib
    INTEGER, INTENT(OUT), OPTIONAL           :: iprint
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: rcut
    LOGICAL, INTENT(OUT), OPTIONAL           :: ldiscrete
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: discrete_step, pmavbmc
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: pbias
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: avbmc_atom
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: avbmc_rmin, avbmc_rmax, &
                                                rmdihedral
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: input_file
    TYPE(mc_molecule_info_type), OPTIONAL, &
      POINTER                                :: mc_molecule_info
    REAL(dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: pmswap_mol
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: pmavbmc_mol, pmtrans_mol, &
                                                pmrot_mol, pmtraion_mol
    TYPE(mc_input_file_type), OPTIONAL, &
      POINTER                                :: mc_input_file, mc_bias_file
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: pmvol_box
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: virial_temps
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: exp_min_val, exp_max_val, &
                                                min_val, max_val
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: eta
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: pmhmc, pmhmc_box
    LOGICAL, INTENT(OUT), OPTIONAL           :: lhmc
    INTEGER, INTENT(OUT), OPTIONAL           :: rand2skip

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_mc_par', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( nstep ) ) nstep = mc_par % nstep
    IF ( PRESENT ( nvirial ) ) nvirial = mc_par % nvirial
    IF ( PRESENT ( iuptrans ) ) iuptrans = mc_par % iuptrans
    IF ( PRESENT ( iupvolume ) ) iupvolume = mc_par % iupvolume
    IF ( PRESENT ( nmoves ) ) nmoves = mc_par % nmoves
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

    IF ( PRESENT ( virial_temps ) ) virial_temps => mc_par % virial_temps
    IF ( PRESENT ( avbmc_atom ) ) avbmc_atom => mc_par % avbmc_atom
    IF ( PRESENT ( avbmc_rmin ) ) avbmc_rmin => mc_par % avbmc_rmin
    IF ( PRESENT ( avbmc_rmax ) ) avbmc_rmax => mc_par % avbmc_rmax
    IF ( PRESENT ( rcut ) ) rcut = mc_par % rcut
    IF ( PRESENT ( discrete_step ) ) discrete_step = mc_par % discrete_step
    IF ( PRESENT ( rmvolume ) ) rmvolume = mc_par % rmvolume
    IF ( PRESENT ( temperature ) ) temperature = mc_par % temperature
    IF ( PRESENT ( pressure ) ) pressure = mc_par % pressure
    IF ( PRESENT ( BETA ) ) BETA = mc_par % BETA
    IF ( PRESENT ( exp_max_val ) ) exp_max_val = mc_par % exp_max_val
    IF ( PRESENT ( exp_min_val ) ) exp_min_val = mc_par % exp_min_val
    IF ( PRESENT ( max_val ) ) max_val = mc_par % max_val
    IF ( PRESENT ( min_val ) ) min_val = mc_par % min_val
    IF ( PRESENT ( pmswap ) ) pmswap = mc_par % pmswap
    IF ( PRESENT ( pmvolume) ) pmvolume = mc_par % pmvolume
    IF ( PRESENT ( lhmc) ) lhmc = mc_par % lhmc
    IF ( PRESENT ( pmhmc) ) pmhmc = mc_par % pmhmc
    IF ( PRESENT ( pmhmc_box) ) pmhmc_box = mc_par % pmhmc_box
    IF ( PRESENT ( pmvol_box) ) pmvol_box = mc_par % pmvol_box
    IF ( PRESENT ( pmtraion ) ) pmtraion = mc_par % pmtraion
    IF ( PRESENT ( pmtrans ) ) pmtrans = mc_par % pmtrans
    IF ( PRESENT ( pmavbmc) ) pmavbmc = mc_par % pmavbmc
    IF ( PRESENT ( pmrot_mol) ) pmrot_mol => mc_par % pmrot_mol
    IF ( PRESENT ( pmtrans_mol) ) pmtrans_mol => mc_par % pmtrans_mol
    IF ( PRESENT ( pmtraion_mol) ) pmtraion_mol => mc_par % pmtraion_mol
    IF ( PRESENT ( pmavbmc_mol) ) pmavbmc_mol => mc_par % pmavbmc_mol
    IF ( PRESENT ( pbias) ) pbias => mc_par % pbias

    IF ( PRESENT ( rmbond ) ) rmbond => mc_par % rmbond
    IF ( PRESENT ( rmangle ) ) rmangle => mc_par % rmangle
    IF ( PRESENT ( rmdihedral ) ) rmdihedral => mc_par % rmdihedral
    IF ( PRESENT ( rmrot ) ) rmrot => mc_par % rmrot
    IF ( PRESENT ( rmtrans ) ) rmtrans => mc_par % rmtrans

    IF ( PRESENT ( eta ) ) eta => mc_par % eta

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

    IF ( PRESENT ( input_file ) ) input_file => mc_par % input_file
    IF ( PRESENT ( mc_molecule_info ) ) mc_molecule_info => mc_par % mc_molecule_info
    IF ( PRESENT ( mc_input_file ) ) mc_input_file => mc_par % mc_input_file
    IF ( PRESENT ( mc_bias_file ) ) mc_bias_file => mc_par % mc_bias_file

    IF ( PRESENT ( pmswap_mol ) ) pmswap_mol => mc_par % pmswap_mol
    IF ( PRESENT ( rand2skip ) ) rand2skip = mc_par % rand2skip
  END SUBROUTINE get_mc_par

! *****************************************************************************
 SUBROUTINE get_mc_molecule_info ( mc_molecule_info, nmol_types,nchain_total,nboxes,&
      names,conf_prob,nchains,nunits,mol_type,nunits_tot,in_box,atom_names,&
      mass)

! *****************************************************************************
!> \brief accesses the private elements of the mc_molecule_info_type
!> \param mc_molecule_info the structure you want the parameters for
!>
!>    Suitable for parallel.
!> \author MJM
    TYPE(mc_molecule_info_type), POINTER     :: mc_molecule_info
    INTEGER, INTENT(OUT), OPTIONAL           :: nmol_types, nchain_total, &
                                                nboxes
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), OPTIONAL, POINTER        :: names
    REAL(dp), DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: conf_prob
    INTEGER, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: nchains
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: nunits, mol_type, nunits_tot, &
                                                in_box
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:, :), OPTIONAL, POINTER     :: atom_names
    REAL(dp), DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: mass

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_mc_molecule_info', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( nboxes ) ) nboxes = mc_molecule_info % nboxes
    IF ( PRESENT ( nchain_total ) ) nchain_total = mc_molecule_info % nchain_total
    IF ( PRESENT ( nmol_types ) ) nmol_types = mc_molecule_info % nmol_types

    IF ( PRESENT ( names ) ) names => mc_molecule_info % names
    IF ( PRESENT ( atom_names ) ) atom_names => mc_molecule_info % atom_names

    IF ( PRESENT ( conf_prob ) ) conf_prob => mc_molecule_info % conf_prob
    IF ( PRESENT ( mass ) ) mass => mc_molecule_info % mass

    IF ( PRESENT ( nchains ) ) nchains => mc_molecule_info % nchains
    IF ( PRESENT ( nunits ) ) nunits => mc_molecule_info % nunits
    IF ( PRESENT ( mol_type ) ) mol_type => mc_molecule_info % mol_type
    IF ( PRESENT ( nunits_tot ) ) nunits_tot => mc_molecule_info % nunits_tot
    IF ( PRESENT ( in_box ) ) in_box => mc_molecule_info % in_box

  END SUBROUTINE get_mc_molecule_info

! *****************************************************************************
SUBROUTINE set_mc_par ( mc_par, rm,cl,&
      diff,nstart,rmvolume,rmbond,rmangle,rmdihedral,rmrot,rmtrans,PROGRAM,&
      nmoves,nswapmoves,lstop,temperature,pressure,iuptrans,iupvolume,&
      pmswap,pmvolume,pmtraion,pmtrans,BETA,rcut,iprint,lbias,nstep,&
      lrestart,ldiscrete,discrete_step,pmavbmc,mc_molecule_info,&
      pmavbmc_mol,pmtrans_mol,pmrot_mol,pmtraion_mol,pmswap_mol,&
      avbmc_rmin,avbmc_rmax,avbmc_atom,pbias,ensemble,pmvol_box,eta,&
      mc_input_file,mc_bias_file,exp_max_val,exp_min_val,min_val,max_val,&
      pmhmc,pmhmc_box,lhmc,ionode,source,group,rand2skip)


! *****************************************************************************
!> \brief changes the private elements of the mc_parameters_type
!> \param mc_par the structure mc parameters you want
!>
!>    Suitable for parallel.
!> \author MJM
    TYPE(mc_simpar_type), POINTER            :: mc_par
    INTEGER, INTENT(IN), OPTIONAL            :: rm, cl, diff, nstart
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: rmvolume
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: rmbond, rmangle, rmdihedral, &
                                                rmrot, rmtrans
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: PROGRAM
    INTEGER, INTENT(IN), OPTIONAL            :: nmoves, nswapmoves
    LOGICAL, INTENT(IN), OPTIONAL            :: lstop
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: temperature, pressure
    INTEGER, INTENT(IN), OPTIONAL            :: iuptrans, iupvolume
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: pmswap, pmvolume, pmtraion, &
                                                pmtrans, BETA, rcut
    INTEGER, INTENT(IN), OPTIONAL            :: iprint
    LOGICAL, INTENT(IN), OPTIONAL            :: lbias
    INTEGER, INTENT(IN), OPTIONAL            :: nstep
    LOGICAL, INTENT(IN), OPTIONAL            :: lrestart, ldiscrete
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: discrete_step, pmavbmc
    TYPE(mc_molecule_info_type), OPTIONAL, &
      POINTER                                :: mc_molecule_info
    REAL(dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: pmavbmc_mol, pmtrans_mol, &
                                                pmrot_mol, pmtraion_mol, &
                                                pmswap_mol, avbmc_rmin, &
                                                avbmc_rmax
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: avbmc_atom
    REAL(dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: pbias
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: ensemble
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: pmvol_box
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: eta
    TYPE(mc_input_file_type), OPTIONAL, &
      POINTER                                :: mc_input_file, mc_bias_file
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: exp_max_val, exp_min_val, &
                                                min_val, max_val, pmhmc, &
                                                pmhmc_box
    LOGICAL, INTENT(IN), OPTIONAL            :: lhmc, ionode
    INTEGER, INTENT(IN), OPTIONAL            :: source, group, rand2skip

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_mc_par', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iparm

! These are the only values that change during the course of the simulation
! or are computed outside of this module

    IF ( PRESENT ( nstep ) ) mc_par % nstep = nstep
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

    IF ( PRESENT ( exp_max_val ) ) mc_par % exp_max_val= exp_max_val
    IF ( PRESENT ( exp_min_val ) ) mc_par % exp_min_val= exp_min_val
    IF ( PRESENT ( max_val ) ) mc_par % max_val = max_val
    IF ( PRESENT ( min_val ) ) mc_par % min_val = min_val
    IF ( PRESENT ( BETA ) ) mc_par % BETA = BETA
    IF ( PRESENT ( temperature ) ) mc_par % temperature = temperature
    IF ( PRESENT ( rcut ) ) mc_par % rcut = rcut
    IF ( PRESENT ( pressure ) ) mc_par % pressure = pressure
    IF ( PRESENT ( pmvolume ) ) mc_par % pmvolume = pmvolume
    IF ( PRESENT ( lhmc ) ) mc_par % lhmc = lhmc
    IF ( PRESENT ( pmhmc ) ) mc_par % pmhmc = pmhmc
    IF ( PRESENT ( pmhmc_box ) ) mc_par % pmhmc_box = pmhmc_box
    IF ( PRESENT ( pmvol_box ) ) mc_par % pmvol_box = pmvol_box
    IF ( PRESENT ( pmswap ) ) mc_par % pmswap = pmswap
    IF ( PRESENT ( pmtrans ) ) mc_par % pmtrans = pmtrans
    IF ( PRESENT ( pmtraion ) ) mc_par % pmtraion = pmtraion
    IF ( PRESENT ( pmavbmc ) ) mc_par % pmavbmc = pmavbmc

    IF ( PRESENT ( discrete_step ) ) mc_par % discrete_step = discrete_step
    IF ( PRESENT ( rmvolume ) ) mc_par % rmvolume = rmvolume

    IF ( PRESENT ( mc_input_file ) ) mc_par % mc_input_file => mc_input_file
    IF ( PRESENT ( mc_bias_file ) ) mc_par % mc_bias_file => mc_bias_file

! cannot just change the pointers, because then we have memory leaks
! and the program crashes if we try to release it (in certain cases)
    IF ( PRESENT ( eta ) ) THEN
       DO iparm=1,SIZE(eta)
          mc_par % eta(iparm) = eta(iparm)
       ENDDO
    ENDIF
    IF ( PRESENT ( rmbond ) ) THEN
       DO iparm=1,SIZE(rmbond)
          mc_par % rmbond(iparm) = rmbond(iparm)
       ENDDO
    ENDIF
    IF ( PRESENT ( rmangle ) ) THEN
       DO iparm=1,SIZE(rmangle)
          mc_par % rmangle(iparm) = rmangle(iparm)
       ENDDO
    ENDIF
    IF ( PRESENT ( rmdihedral ) ) THEN
       DO iparm=1,SIZE(rmdihedral)
          mc_par % rmdihedral(iparm) = rmdihedral(iparm)
       ENDDO
    ENDIF
    IF ( PRESENT ( rmrot ) ) THEN
       DO iparm=1,SIZE(rmrot)
          mc_par % rmrot(iparm) = rmrot(iparm)
       ENDDO
    ENDIF
    IF ( PRESENT ( rmtrans ) ) THEN
       DO iparm=1,SIZE(rmtrans)
          mc_par % rmtrans(iparm) = rmtrans(iparm)
       ENDDO
    ENDIF

    IF ( PRESENT ( pmavbmc_mol ) ) THEN
       DO iparm=1,SIZE(pmavbmc_mol)
          mc_par % pmavbmc_mol(iparm) = pmavbmc_mol(iparm)
       ENDDO
    ENDIF
    IF ( PRESENT ( pmtrans_mol ) ) THEN
       DO iparm=1,SIZE(pmtrans_mol)
          mc_par % pmtrans_mol(iparm) = pmtrans_mol(iparm)
       ENDDO
    ENDIF
    IF ( PRESENT ( pmrot_mol ) ) THEN
       DO iparm=1,SIZE(pmrot_mol)
          mc_par % pmrot_mol(iparm) = pmrot_mol(iparm)
       ENDDO
    ENDIF
    IF ( PRESENT ( pmtraion_mol ) ) THEN
       DO iparm=1,SIZE(pmtraion_mol)
          mc_par % pmtraion_mol(iparm) = pmtraion_mol(iparm)
       ENDDO
    ENDIF
    IF ( PRESENT ( pmswap_mol ) ) THEN
       DO iparm=1,SIZE(pmswap_mol)
          mc_par % pmswap_mol(iparm) = pmswap_mol(iparm)
       ENDDO
    ENDIF

    IF ( PRESENT ( avbmc_atom ) ) THEN
       DO iparm=1,SIZE(avbmc_atom)
          mc_par % avbmc_atom(iparm) = avbmc_atom(iparm)
       ENDDO
    ENDIF
    IF ( PRESENT ( avbmc_rmin ) ) THEN
       DO iparm=1,SIZE(avbmc_rmin)
          mc_par % avbmc_rmin(iparm) = avbmc_rmin(iparm)
       ENDDO
    ENDIF
    IF ( PRESENT ( avbmc_rmax ) ) THEN
       DO iparm=1,SIZE(avbmc_rmax)
          mc_par % avbmc_rmax(iparm) = avbmc_rmax(iparm)
       ENDDO
    ENDIF
    IF ( PRESENT ( pbias ) ) THEN
       DO iparm=1,SIZE(pbias)
          mc_par % pbias(iparm) = pbias(iparm)
       ENDDO
    ENDIF

    IF ( PRESENT ( PROGRAM ) ) mc_par % program = PROGRAM
    IF ( PRESENT ( ensemble ) ) mc_par % ensemble = ensemble

    IF(PRESENT(mc_molecule_info)) mc_par%mc_molecule_info => mc_molecule_info
    IF(PRESENT(source)) mc_par%source=source
    IF(PRESENT(group)) mc_par%group=group
    IF(PRESENT(ionode)) mc_par%ionode=ionode

    IF(PRESENT(rand2skip)) mc_par%rand2skip=rand2skip
  END SUBROUTINE set_mc_par

! *****************************************************************************
!> \brief creates (allocates) the mc_simulation_parameters type
!> \param mc_par the structure that will be created (allocated)
!> \param nmol_types the number of molecule types in the system
!> \author MJM
! *****************************************************************************
  SUBROUTINE mc_sim_par_create ( mc_par ,nmol_types)

    TYPE(mc_simpar_type), POINTER            :: mc_par
    INTEGER, INTENT(IN)                      :: nmol_types

    CHARACTER(len=*), PARAMETER :: routineN = 'mc_sim_par_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    NULLIFY(mc_par)

    ALLOCATE(mc_par)
    ALLOCATE (mc_par%pmtraion_mol(1:nmol_types),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%pmtraion_mol",nmol_types*dp_size)
    ALLOCATE (mc_par%pmtrans_mol(1:nmol_types),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%pmtrans_mol",nmol_types*dp_size)
    ALLOCATE (mc_par%pmrot_mol(1:nmol_types),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%pmrot_mol",nmol_types*dp_size)
    ALLOCATE (mc_par%pmswap_mol(1:nmol_types),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%pmswap_mol",nmol_types*dp_size)

    ALLOCATE (mc_par%eta(1:nmol_types),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%eta",nmol_types*dp_size)

    ALLOCATE (mc_par%rmbond(1:nmol_types),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%rmbond",nmol_types*dp_size)
    ALLOCATE (mc_par%rmangle(1:nmol_types),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%rmangle",nmol_types*dp_size)
    ALLOCATE (mc_par%rmdihedral(1:nmol_types),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%rmdihedral",nmol_types*dp_size)
    ALLOCATE (mc_par%rmtrans(1:nmol_types),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%rmtrans",nmol_types*dp_size)
    ALLOCATE (mc_par%rmrot(1:nmol_types),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%rmrot",nmol_types*dp_size)

    ALLOCATE (mc_par%avbmc_atom(1:nmol_types),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%avbmc_atom",nmol_types*int_size)
    ALLOCATE (mc_par%avbmc_rmin(1:nmol_types),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%avbmc_rmin",nmol_types*dp_size)
    ALLOCATE (mc_par%avbmc_rmax(1:nmol_types),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%avbmc_rmax",nmol_types*dp_size)
    ALLOCATE (mc_par%pmavbmc_mol(1:nmol_types),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%pmavbmc_mol",nmol_types*dp_size)
    ALLOCATE (mc_par%pbias(1:nmol_types),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%pbias",nmol_types*dp_size)

    ALLOCATE (mc_par%mc_input_file,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%mc_input_file",1)
    ALLOCATE (mc_par%mc_bias_file,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%mc_bias_file",1)

  END SUBROUTINE mc_sim_par_create

! *****************************************************************************
!> \brief destroys (deallocates) the mc_simulation_parameters type
!> \param mc_par the structure that will be destroyed
!> \author MJM
! *****************************************************************************
  SUBROUTINE mc_sim_par_destroy ( mc_par )

    TYPE(mc_simpar_type), POINTER            :: mc_par

    CHARACTER(len=*), PARAMETER :: routineN = 'mc_sim_par_destroy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    DEALLOCATE (mc_par%mc_input_file,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%mc_input_file")
    DEALLOCATE (mc_par%mc_bias_file,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%mc_bias_file")

   DEALLOCATE (mc_par%pmtraion_mol,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%pmtraion_mol")
    DEALLOCATE (mc_par%pmtrans_mol,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%pmtrans_mol")
    DEALLOCATE (mc_par%pmrot_mol,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%pmrot_mol")
    DEALLOCATE (mc_par%pmswap_mol,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%pmswap_mol")

    DEALLOCATE (mc_par%eta,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%eta")

    DEALLOCATE (mc_par%rmbond,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%rmbond")
    DEALLOCATE (mc_par%rmangle,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%rmangle")
    DEALLOCATE (mc_par%rmdihedral,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%rmdihedral")
    DEALLOCATE (mc_par%rmtrans,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%rmtrans")
    DEALLOCATE (mc_par%rmrot,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%rmrot")

    DEALLOCATE (mc_par%avbmc_atom,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%avbmc_atom")
    DEALLOCATE (mc_par%avbmc_rmin,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%avbmc_rmin")
    DEALLOCATE (mc_par%avbmc_rmax,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%avbmc_rmax")
    DEALLOCATE (mc_par%pmavbmc_mol,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%pmavbmc_mol")
    DEALLOCATE (mc_par%pbias,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%pbias")
   IF(mc_par%ensemble == "VIRIAL") THEN
    DEALLOCATE (mc_par%virial_temps,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par%virial_temps")
   END IF

   DEALLOCATE (mc_par,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par")
    NULLIFY(mc_par)

 END SUBROUTINE mc_sim_par_destroy

! *****************************************************************************
!> \brief creates (allocates) the mc_input_file_type
!> \param mc_input_file the structure that will be created (allocated)
!> \param input_file_name the name of the file to read
!> \param nmol_types the number of molecule types in the system
!> \author MJM
! *****************************************************************************
  SUBROUTINE mc_input_file_create (mc_input_file,input_file_name,&
       mc_molecule_info,empty_coords,lhmc,error)

    TYPE(mc_input_file_type), POINTER        :: mc_input_file
    CHARACTER(LEN=*), INTENT(IN)             :: input_file_name
    TYPE(mc_molecule_info_type), POINTER     :: mc_molecule_info
    REAL(dp), DIMENSION(:, :)                :: empty_coords
    LOGICAL, INTENT(IN)                      :: lhmc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mc_input_file_create', &
      routineP = moduleN//':'//routineN

    CHARACTER(default_string_length)         :: cdum, line, method_name
    CHARACTER(default_string_length), &
      DIMENSION(:, :), POINTER               :: atom_names
    INTEGER :: abc_column, abc_row, cell_column, cell_row, idum, iline, &
      io_stat, irep, istat, iunit, iw, line_number, nlines, nmol_types, &
      nstart, row_number, unit
    INTEGER, DIMENSION(:), POINTER           :: nunits
    TYPE(cp_logger_type), POINTER            :: logger

! some stuff in case we need to write error messages to the screen

    logger=>cp_error_get_logger(error)
    iw = cp_logger_get_default_io_unit(logger)

! allocate the array we'll need in case we have an empty box
    CALL get_mc_molecule_info(mc_molecule_info,nmol_types=nmol_types,&
         nunits=nunits,atom_names=atom_names)
    ALLOCATE (mc_input_file%coordinates_empty(1:3,1:nunits(1)),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_input_file%coordinates_empty",nunits(1)*5)
    ALLOCATE (mc_input_file%atom_names_empty(1:nunits(1)),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_input_file%atom_names_empty",nunits(1)*5)
    DO iunit=1,nunits(1)
       mc_input_file%atom_names_empty(iunit)=atom_names(iunit,1)
       mc_input_file%coordinates_empty(1:3,iunit)=empty_coords(1:3,iunit)
    ENDDO
    mc_input_file%nunits_empty=nunits(1)

! allocate some arrays
    ALLOCATE (mc_input_file%mol_set_nmol_row(1:nmol_types),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_input_file%mol_set_nmol_row",nmol_types*int_size)
    ALLOCATE (mc_input_file%mol_set_nmol_column(1:nmol_types),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_input_file%mol_set_nmol_column",nmol_types*int_size)

! now read in and store the input file, first finding out how many lines it is
    CALL open_file(file_name=input_file_name,&
         unit_number=unit,file_action='READ',file_status='OLD')

    nlines=0
    DO
       READ(unit,'(A)',IOSTAT=io_stat) line
       IF(io_stat .NE. 0) EXIT
       nlines=nlines+1
    ENDDO

    ALLOCATE (mc_input_file%text(1:nlines),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_input_file%text",nlines*default_string_length)

    REWIND(unit)
    DO iline=1,nlines
       READ(unit,'(A)') mc_input_file%text(iline)
    ENDDO

! close the input file
    CALL close_file(unit_number=unit)

! now we need to find the row and column numbers of a variety of things
! first, Let's find the first END after GLOBAL
    CALL mc_parse_text(mc_input_file%text,1,nlines,"&GLOBAL",.TRUE.,&
         mc_input_file%global_row_end,idum)
    IF(mc_input_file%global_row_end == 0) THEN
       IF (iw>0) THEN
          WRITE(iw,*)
          WRITE(iw,*) 'File name ',input_file_name
       END IF
       CALL stop_program(routineN,moduleN,__LINE__,&
            'Could not find &END after &GLOBAL (make sure & is in the same column) and last line is blank')
    ENDIF

! Let's find the first END after MOTION...we might need this whole section
! because of hybrid MD/MC moves, which requires the MD information to always
! be attached to the force env
    CALL mc_parse_text(mc_input_file%text,1,nlines,"&MOTION",.TRUE.,&
         mc_input_file%motion_row_end,idum,start_row_number=mc_input_file%motion_row_start)
    IF(mc_input_file%motion_row_end == 0) THEN
       IF (iw>0) THEN
          WRITE(iw,*)
          WRITE(iw,*) 'File name ',input_file_name,mc_input_file%motion_row_start
       END IF
       CALL stop_program(routineN,moduleN,__LINE__,&
            'Could not find &END after &MOTION (make sure & is in the same column) and last line is blank')
    ENDIF

! first, let's find the first END after &COORD, and the line of &COORD
    CALL mc_parse_text(mc_input_file%text,1,nlines,"&COORD",.TRUE.,&
         mc_input_file%coord_row_end,idum)
    CALL mc_parse_text(mc_input_file%text,1,nlines,"&COORD",.FALSE.,&
         mc_input_file%coord_row_start,idum)

    IF(mc_input_file%coord_row_end == 0) THEN
       IF (iw>0) THEN
          WRITE(iw,*)
          WRITE(iw,*) 'File name ',input_file_name
       END IF
       CALL stop_program(routineN,moduleN,__LINE__,&
            'Could not find &END after &COORD (make sure & is the first in the same column after &COORD)')
    ENDIF

! now the RUN_TYPE
    CALL mc_parse_text(mc_input_file%text,1,nlines,"RUN_TYPE",.FALSE.,&
         mc_input_file%run_type_row, mc_input_file%run_type_column)
    mc_input_file%run_type_column=mc_input_file%run_type_column+9
    IF(mc_input_file%run_type_row == 0) THEN
       IF (iw>0) THEN
          WRITE(iw,*)
          WRITE(iw,*) 'File name ',input_file_name
       END IF
       CALL stop_program(routineN,moduleN,__LINE__,&
            'Could not find RUN_TYPE in the input file (should be in the &GLOBAL section).')
    ENDIF

! now the CELL stuff..we don't care about REF_CELL since that should
! never change if it's there
    CALL mc_parse_text(mc_input_file%text,1,nlines,"&CELL",.FALSE.,&
         mc_input_file%cell_row, mc_input_file%cell_column)
! now find the ABC input line after CELL
    CALL mc_parse_text(mc_input_file%text,mc_input_file%cell_row+1,nlines,&
         "ABC",.FALSE.,abc_row,abc_column)
! is there a &CELL inbetween?  If so, that ABC will be for the ref_cell
! and we need to find the next one
    CALL mc_parse_text(mc_input_file%text,mc_input_file%cell_row+1,abc_row,&
         "&CELL",.FALSE.,cell_row,cell_column)
    IF(cell_row == 0) THEN ! nothing in between...we found the correct ABC
       mc_input_file%cell_row=abc_row
       mc_input_file%cell_column=abc_column+4
    ELSE
       CALL mc_parse_text(mc_input_file%text,abc_row+1,nlines,&
            "ABC",.FALSE.,mc_input_file%cell_row, mc_input_file%cell_column)
    ENDIF
    IF(mc_input_file%cell_row == 0) THEN
       IF (iw>0) THEN
          WRITE(iw,*)
          WRITE(iw,*) 'File name ',input_file_name
       END IF
       CALL stop_program(routineN,moduleN,__LINE__,&
            'Could not find &CELL section in the input file.  Or could not find the ABC line after it.')
    ENDIF

! now we need all the MOL_SET NMOLS indicies
    IF(.NOT. lhmc)THEN
       irep=0
       nstart=1
       DO
          CALL mc_parse_text(mc_input_file%text,nstart,nlines,"&MOLECULE",&
               .FALSE.,row_number, idum)
          IF(row_number == 0) EXIT
          nstart=row_number+1
          irep=irep+1
          CALL mc_parse_text(mc_input_file%text,nstart,nlines,"NMOL",&
               .FALSE.,mc_input_file%mol_set_nmol_row(irep), &
               mc_input_file%mol_set_nmol_column(irep))
          mc_input_file%mol_set_nmol_column(irep)=mc_input_file%mol_set_nmol_column(irep)+5

       ENDDO
       IF(irep .NE. nmol_types) THEN
          IF (iw>0) THEN
             WRITE(iw,*)
             WRITE(iw,*) 'File name ',input_file_name
             WRITE(iw,*) 'Number of molecule types ',nmol_types
             WRITE(iw,*) '&MOLECULE sections found ',irep
          END IF
          CALL stop_program(routineN,moduleN,__LINE__,&
               'Did not find MOLECULE sections for every molecule in the simulation (make sure both input files have all types)')
       ENDIF
    ENDIF

! now let's find the type of force_env...right now, I'm only concerned with
! QS, and Fist, though it should be trivial to add others
    CALL mc_parse_text(mc_input_file%text,1,nlines,&
         "METHOD",.FALSE.,line_number, idum)
    READ(mc_input_file%text(line_number),*) cdum,method_name
    CALL uppercase(method_name)
    SELECT CASE(TRIM(ADJUSTL(method_name)))
    CASE("FIST")
       mc_input_file%in_use=use_fist_force
    CASE("QS","QUICKSTEP")
       mc_input_file%in_use=use_qs_force
    CASE default
       CALL stop_program(routineN,moduleN,__LINE__,&
            'Cannot determine the type of force_env we are using (check METHOD)')
    END SELECT

  END SUBROUTINE mc_input_file_create

! *****************************************************************************
!> \brief destroys (deallocates) things in the mc_input_file_type
!> \param mc_input_file the structure that will be released (deallocated)
!> \author MJM
! *****************************************************************************
  SUBROUTINE mc_input_file_destroy(mc_input_file)

    TYPE(mc_input_file_type), POINTER        :: mc_input_file

    CHARACTER(len=*), PARAMETER :: routineN = 'mc_input_file_destroy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    DEALLOCATE (mc_input_file%mol_set_nmol_row,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_input_file%mol_set_nmol_row")
    DEALLOCATE (mc_input_file%mol_set_nmol_column,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_input_file%mol_set_nmol_column")
    DEALLOCATE (mc_input_file%text,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_input_file%text")
    DEALLOCATE (mc_input_file%atom_names_empty,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_input_file%atom_names_empty")
    DEALLOCATE (mc_input_file%coordinates_empty,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_input_file%coordinates_empty")

  END SUBROUTINE mc_input_file_destroy

! *****************************************************************************
!> \brief a basic text parser used to find the row and column numbers of various strings
!>      in the input file, to store as indices for when we create a dat_file...
!>      returns 0 for the row if nothing is found
!> \param text the text to parse
!> \param nstart the line to start searching from
!> \param nend the line to end searching
!> \param string_search the text we're looking for
!> \param lend if .TRUE., find the &END that comes after string_search...assumes that
!>            the row is the first row with & in the same column as the search string
!> \param row_number the row the string is first found on
!> \param column_number the column number that the string first appears on
!> \author MJM
! *****************************************************************************
SUBROUTINE mc_parse_text(text,nstart,nend,string_search,lend,&
     row_number,column_number,start_row_number)

    CHARACTER(LEN=*), DIMENSION(:), &
      INTENT(IN)                             :: text
    INTEGER, INTENT(IN)                      :: nstart, nend
    CHARACTER(LEN=*), INTENT(IN)             :: string_search
    LOGICAL, INTENT(IN)                      :: lend
    INTEGER, INTENT(OUT)                     :: row_number, column_number
    INTEGER, INTENT(OUT), OPTIONAL           :: start_row_number

    CHARACTER(default_string_length)         :: text_temp
    INTEGER                                  :: iline, index_string, jline, &
                                                string_size

! did not see any string utilities in the code to help, so here I go

  row_number=0
  column_number=0

! how long is our string to search for?
  string_size=LEN_TRIM(string_search)
  whole_file:DO iline=nstart,nend

     index_string=-1
     index_string=INDEX(text(iline),string_search(1:string_size))

     IF(index_string .GT. 0) THEN ! we found one
        IF(.NOT. lend) THEN
           row_number=iline
           column_number=index_string
        ELSE
           IF(PRESENT(start_row_number)) start_row_number=iline
           column_number=index_string
           DO jline=iline+1,nend
! now we find the &END that matches up with this one...
! I need proper indentation because I'm not very smart
              text_temp=text(jline)
              IF(text_temp(column_number:column_number) .EQ. "&") THEN
                 row_number=jline
                 EXIT
              ENDIF
           ENDDO
        ENDIF
        EXIT whole_file
     ENDIF
  ENDDO whole_file

END SUBROUTINE mc_parse_text

! *****************************************************************************
!> \brief reads in the Monte Carlo simulation parameters from an input file
!> \param mc_par the structure that will store the parameters
!> \param globenv the global environment for the simulation
!> \param input_file_name the name of the input_file
!> \param input_file the structure that contains all the keywords in the input file
!> \param force_env_section used to grab the type of force_env
!> \param error the cp_error_type
!> \author MJM
! *****************************************************************************
SUBROUTINE read_mc_section ( mc_par, para_env, globenv, input_file_name, input_file, force_env_section,&
     error )

    TYPE(mc_simpar_type), POINTER            :: mc_par
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    CHARACTER(LEN=*), INTENT(IN)             :: input_file_name
    TYPE(section_vals_type), POINTER         :: input_file, force_env_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_mc_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=5)                         :: box_string, molecule_string, &
                                                tab_box_string, tab_string, &
                                                tab_string_int
    CHARACTER(LEN=default_string_length)     :: format_box_string, &
                                                format_string, &
                                                format_string_int
    INTEGER                                  :: handle, ia, ie, imol, istat, &
                                                itype, iw, method_name_id, &
                                                nmol_types, stop_num
    INTEGER, DIMENSION(:), POINTER           :: temp_i_array
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:), POINTER          :: temp_r_array
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: mc_section

! begin the timing of the subroutine

    CPPrecondition(ASSOCIATED(input_file),cp_failure_level,routineP,error,failure)

      CALL timeset(routineN,handle)

      NULLIFY(mc_section)
      mc_section => section_vals_get_subs_vals(input_file,&
         "MOTION%MC",error=error)

! need the input file sturcutre that we're reading from for when we make
! dat files
      mc_par % input_file => input_file

! set the ionode and mepos
      mc_par % ionode = para_env%ionode
      mc_par % group = para_env%group
      mc_par % source = para_env%source

! for convenience
      nmol_types=mc_par%mc_molecule_info%nmol_types

!..defaults...most are set in input_cp2k_motion
      mc_par % nstart = 0
      CALL section_vals_val_get(force_env_section,"METHOD",i_val=method_name_id,error=error)
      SELECT CASE (method_name_id)
      CASE (do_fist)
         mc_par % iprint = 100
      CASE (do_qs)
         mc_par % iprint = 1
      END SELECT

      logger=>cp_error_get_logger(error)
      iw = cp_logger_get_default_io_unit(logger)
      IF(iw>0) WRITE ( iw, * )

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
      mc_par % fft_lib = globenv%default_fft_library

! get all the move probabilities first...if we are not doing certain types of moves, we don't care
! if the values for those move parameters are strange

! find out if we're only doing HMC moves...we can ignore a lot of extra information
! then, which would ordinarily be cause for concern
      CALL section_vals_val_get(mc_section,"MOVE_PROBABILITIES%PMHMC",r_val=mc_par%pmhmc,error=error)
      CALL section_vals_val_get(mc_section,"MOVE_PROBABILITIES%PMSWAP",r_val=mc_par%pmswap,error=error)
      CALL section_vals_val_get(mc_section,"MOVE_PROBABILITIES%PMVOLUME",r_val=mc_par%pmvolume,error=error)
      CALL section_vals_val_get(mc_section,"MOVE_PROBABILITIES%PMAVBMC",r_val=mc_par%pmavbmc,error=error)
      CALL section_vals_val_get(mc_section,"MOVE_PROBABILITIES%PMTRANS",r_val=mc_par%pmtrans,error=error)
      CALL section_vals_val_get(mc_section,"MOVE_PROBABILITIES%PMTRAION",r_val=mc_par%pmtraion,error=error)


 ! first, grab all the integer values
      CALL section_vals_val_get(mc_section,"NSTEP",i_val=mc_par%nstep,error=error)
      CALL section_vals_val_get(mc_section,"NMOVES",i_val=mc_par%nmoves,error=error)
      CALL section_vals_val_get(mc_section,"NSWAPMOVES",i_val=mc_par%nswapmoves,error=error)
      CALL section_vals_val_get(mc_section,"MOVE_UPDATES%IUPVOLUME",&
           i_val=mc_par%iupvolume,error=error)
      CALL section_vals_val_get(mc_section,"NVIRIAL",i_val=mc_par%nvirial,error=error)
      CALL section_vals_val_get(mc_section,"MOVE_UPDATES%IUPTRANS",&
           i_val=mc_par%iuptrans,error=error)
      CALL section_vals_val_get(mc_section,"IPRINT",i_val=mc_par%iprint,error=error)
! now an integer array
      CALL section_vals_val_get(mc_section,"AVBMC%AVBMC_ATOM",i_vals=temp_i_array,error=error)


      IF(mc_par%pmhmc - mc_par%pmswap >= 1.0_dp .AND. mc_par%pmhmc == 1.0_dp)THEN
         mc_par%lhmc=.TRUE.
      ELSE
         mc_par%lhmc=.FALSE.
      ENDIF

      IF( .NOT. mc_par%lhmc)THEN
         IF(SIZE(temp_i_array) .NE. nmol_types)THEN
            CALL stop_program(routineN,moduleN,__LINE__,&
                 'AVBMC_ATOM must have as many values as there are molecule types.')
         ENDIF
      ENDIF
      DO imol=1,SIZE(temp_i_array)
         mc_par%avbmc_atom(imol)=temp_i_array(imol)
      ENDDO

! now the real values
      CALL section_vals_val_get(mc_section,"PRESSURE",r_val=mc_par%pressure,error=error)
      CALL section_vals_val_get(mc_section,"TEMPERATURE",r_val=mc_par%temperature,error=error)
      CALL section_vals_val_get(mc_section,"MOVE_PROBABILITIES%BOX_PROBABILITIES%PMHMC_BOX",&
           r_val=mc_par%pmhmc_box,error=error)
      CALL section_vals_val_get(mc_section,"MOVE_PROBABILITIES%BOX_PROBABILITIES%PMVOL_BOX",&
           r_val=mc_par%pmvol_box,error=error)
      CALL section_vals_val_get(mc_section,"DISCRETE_STEP",r_val=mc_par%discrete_step,error=error)
      CALL section_vals_val_get(mc_section,"MAX_DISPLACEMENTS%BOX_DISPLACEMENTS%RMVOLUME",&
           r_val=mc_par%rmvolume,error=error)

! finally the character values
      CALL section_vals_val_get(mc_section,"ENSEMBLE",c_val=mc_par%ensemble,error=error)
      CALL section_vals_val_get(mc_section,"RESTART_FILE_NAME",c_val=mc_par%restart_file_name,error=error)
      CALL section_vals_val_get(mc_section,"COORDINATE_FILE_NAME",c_val=mc_par%coords_file,error=error)
      CALL section_vals_val_get(mc_section,"ENERGY_FILE_NAME",c_val=mc_par%energy_file,error=error)
      CALL section_vals_val_get(mc_section,"MOVES_FILE_NAME",c_val=mc_par%moves_file,error=error)
      CALL section_vals_val_get(mc_section,"MOLECULES_FILE_NAME",c_val=mc_par%molecules_file,error=error)
      CALL section_vals_val_get(mc_section,"CELL_FILE_NAME",c_val=mc_par%cell_file,error=error)
      CALL section_vals_val_get(mc_section,"DATA_FILE_NAME",c_val=mc_par%data_file,error=error)
      CALL section_vals_val_get(mc_section,"MAX_DISP_FILE_NAME",c_val=mc_par%displacement_file,error=error)
      CALL section_vals_val_get(mc_section,"BOX2_FILE_NAME",c_val=mc_par%box2_file,error=error)

! set the values of the arrays...if we just point, we have problems when we start fooling around
! with releasing force_envs and wonky values get set (despite that these are private)
      IF(mc_par%ensemble == "VIRIAL") THEN
        CALL section_vals_val_get(mc_section,"VIRIAL_TEMPS",r_vals=temp_r_array,error=error)
! yes, I'm allocating here...I cannot find a better place to do it, though
        ALLOCATE (mc_par%virial_temps(1:SIZE(temp_r_array)),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
           "virial_temps",dp_size*SIZE(temp_r_array))
        DO imol=1,SIZE(temp_r_array)
           mc_par%virial_temps(imol)=temp_r_array(imol)
        ENDDO
      END IF
! all of these arrays should have one value for each type of molecule...so check that
      CALL section_vals_val_get(mc_section,"AVBMC%AVBMC_RMIN",r_vals=temp_r_array,error=error)
      IF( .NOT. mc_par%lhmc)THEN
         IF(SIZE(temp_r_array) .NE. nmol_types)THEN
            CALL stop_program(routineN,moduleN,__LINE__,&
                 'AVBMC_RMIN must have as many values as there are molecule types.')
         ENDIF
      ENDIF
      DO imol=1,SIZE(temp_r_array)
         mc_par%avbmc_rmin(imol)=temp_r_array(imol)
      ENDDO
      CALL section_vals_val_get(mc_section,"AVBMC%AVBMC_RMAX",r_vals=temp_r_array,error=error)
      IF( .NOT. mc_par%lhmc)THEN
         IF(SIZE(temp_r_array) .NE. nmol_types)THEN
            CALL stop_program(routineN,moduleN,__LINE__,&
                 'AVBMC_RMAXL must have as many values as there are molecule types.')
         ENDIF
      ENDIF
      DO imol=1,SIZE(temp_r_array)
         mc_par%avbmc_rmax(imol)=temp_r_array(imol)
      ENDDO
      CALL section_vals_val_get(mc_section,"AVBMC%PBIAS",r_vals=temp_r_array,error=error)
      IF( .NOT. mc_par%lhmc)THEN
         IF(SIZE(temp_r_array) .NE. nmol_types)THEN
            CALL stop_program(routineN,moduleN,__LINE__,&
                 'PBIAS must have as many values as there are molecule types.')
         ENDIF
      ENDIF
      DO imol=1,SIZE(temp_r_array)
         mc_par%pbias(imol)=temp_r_array(imol)
      ENDDO
      CALL section_vals_val_get(mc_section,"MOVE_PROBABILITIES%MOL_PROBABILITIES%PMAVBMC_MOL",&
           r_vals=temp_r_array,error=error)
      IF( .NOT. mc_par%lhmc)THEN
         IF(SIZE(temp_r_array) .NE. nmol_types)THEN
            CALL stop_program(routineN,moduleN,__LINE__,&
                 'PMAVBMC_MOL must have as many values as there are molecule types.')
         ENDIF
      ENDIF
      DO imol=1,SIZE(temp_r_array)
         mc_par%pmavbmc_mol(imol)=temp_r_array(imol)
      ENDDO
      CALL section_vals_val_get(mc_section,"ETA",r_vals=temp_r_array,error=error)
      IF( .NOT. mc_par%lhmc)THEN
         IF(SIZE(temp_r_array) .NE. nmol_types)THEN
            CALL stop_program(routineN,moduleN,__LINE__,&
                 'ETA must have as many values as there are molecule types.')
         ENDIF
      ENDIF
      DO imol=1,SIZE(temp_r_array)
         mc_par%eta(imol)=temp_r_array(imol)
      ENDDO
      CALL section_vals_val_get(mc_section,"MAX_DISPLACEMENTS%MOL_DISPLACEMENTS%RMBOND",&
           r_vals=temp_r_array,error=error)
      IF( .NOT. mc_par%lhmc)THEN
         IF(SIZE(temp_r_array) .NE. nmol_types)THEN
            CALL stop_program(routineN,moduleN,__LINE__,&
                 'RMBOND must have as many values as there are molecule types.')
         ENDIF
      ENDIF
      DO imol=1,SIZE(temp_r_array)
         mc_par%rmbond(imol)=temp_r_array(imol)
      ENDDO
      CALL section_vals_val_get(mc_section,"MAX_DISPLACEMENTS%MOL_DISPLACEMENTS%RMANGLE",&
           r_vals=temp_r_array,error=error)
      IF( .NOT. mc_par%lhmc)THEN
         IF(SIZE(temp_r_array) .NE. nmol_types)THEN
            CALL stop_program(routineN,moduleN,__LINE__,&
                 'RMANGLE must have as many values as there are molecule types.')
         ENDIF
      ENDIF
      DO imol=1,SIZE(temp_r_array)
         mc_par%rmangle(imol)=temp_r_array(imol)
      ENDDO
      CALL section_vals_val_get(mc_section,"MAX_DISPLACEMENTS%MOL_DISPLACEMENTS%RMDIHEDRAL",&
           r_vals=temp_r_array,error=error)
      IF( .NOT. mc_par%lhmc)THEN
         IF(SIZE(temp_r_array) .NE. nmol_types)THEN
            CALL stop_program(routineN,moduleN,__LINE__,&
                 'RMDIHEDRAL must have as many values as there are molecule types.')
         ENDIF
      ENDIF
      DO imol=1,SIZE(temp_r_array)
         mc_par%rmdihedral(imol)=temp_r_array(imol)
      ENDDO
      CALL section_vals_val_get(mc_section,"MAX_DISPLACEMENTS%MOL_DISPLACEMENTS%RMROT",&
           r_vals=temp_r_array,error=error)
      IF( .NOT. mc_par%lhmc)THEN
         IF(SIZE(temp_r_array) .NE. nmol_types)THEN
            CALL stop_program(routineN,moduleN,__LINE__,&
                 'RMROT must have as many values as there are molecule types.')
         ENDIF
      ENDIF
      DO imol=1,SIZE(temp_r_array)
         mc_par%rmrot(imol)=temp_r_array(imol)
      ENDDO
      CALL section_vals_val_get(mc_section,"MAX_DISPLACEMENTS%MOL_DISPLACEMENTS%RMTRANS",&
           r_vals=temp_r_array,error=error)
      IF( .NOT. mc_par%lhmc)THEN
         IF(SIZE(temp_r_array) .NE. nmol_types)THEN
            CALL stop_program(routineN,moduleN,__LINE__,&
                 'RMTRANS must have as many values as there are molecule types.')
         ENDIF
      ENDIF
      DO imol=1,SIZE(temp_r_array)
         mc_par%rmtrans(imol)=temp_r_array(imol)
      ENDDO

      CALL section_vals_val_get(mc_section,"MOVE_PROBABILITIES%MOL_PROBABILITIES%PMTRAION_MOL",&
           r_vals=temp_r_array,error=error)
      IF( .NOT. mc_par%lhmc)THEN
         IF(SIZE(temp_r_array) .NE. nmol_types)THEN
            CALL stop_program(routineN,moduleN,__LINE__,&
                 'PMTRAION_MOL must have as many values as there are molecule types.')
         ENDIF
      ENDIF
      DO imol=1,SIZE(temp_r_array)
         mc_par%pmtraion_mol(imol)=temp_r_array(imol)
      ENDDO

      CALL section_vals_val_get(mc_section,"MOVE_PROBABILITIES%MOL_PROBABILITIES%PMTRANS_MOL",&
           r_vals=temp_r_array,error=error)
      IF( .NOT. mc_par%lhmc)THEN
         IF(SIZE(temp_r_array) .NE. nmol_types)THEN
            CALL stop_program(routineN,moduleN,__LINE__,&
                 'PMTRANS_MOL must have as many values as there are molecule types.')
         ENDIF
      ENDIF
      DO imol=1,SIZE(temp_r_array)
         mc_par%pmtrans_mol(imol)=temp_r_array(imol)
      ENDDO

      CALL section_vals_val_get(mc_section,"MOVE_PROBABILITIES%MOL_PROBABILITIES%PMROT_MOL",&
           r_vals=temp_r_array,error=error)
      IF( .NOT. mc_par%lhmc)THEN
         IF(SIZE(temp_r_array) .NE. nmol_types)THEN
            CALL stop_program(routineN,moduleN,__LINE__,&
                 'PMROT_MOL must have as many values as there are molecule types.')
         ENDIF
      ENDIF
      DO imol=1,SIZE(temp_r_array)
         mc_par%pmrot_mol(imol)=temp_r_array(imol)
      ENDDO

      CALL section_vals_val_get(mc_section,"MOVE_PROBABILITIES%MOL_PROBABILITIES%PMSWAP_MOL",&
           r_vals=temp_r_array,error=error)
      IF( .NOT. mc_par%lhmc)THEN
         IF(SIZE(temp_r_array) .NE. nmol_types)THEN
            CALL stop_program(routineN,moduleN,__LINE__,&
                 'PMSWAP_MOL must have as many values as there are molecule types.')
         ENDIF
      ENDIF
      DO imol=1,SIZE(temp_r_array)
         mc_par%pmswap_mol(imol)=temp_r_array(imol)
      ENDDO

! now some logical values
      CALL section_vals_val_get(mc_section,"LBIAS",l_val=mc_par%lbias,error=error)
      CALL section_vals_val_get(mc_section,"LDISCRETE",l_val=mc_par%ldiscrete,error=error)
      CALL section_vals_val_get(mc_section,"LSTOP",l_val=mc_par%lstop,error=error)
      CALL section_vals_val_get(mc_section,"RESTART",l_val=mc_par%lrestart,error=error)
      CALL section_vals_val_get(mc_section,"LBIAS",l_val=mc_par%lbias,error=error)

!..end of parsing the input section

!..write some information to output
      IF (iw>0) THEN
         WRITE(box_string,'(I2)') mc_par%mc_molecule_info%nboxes
         WRITE(molecule_string,'(I2)') mc_par%mc_molecule_info%nmol_types
         WRITE(tab_string,'(I4)') 81-10*mc_par%mc_molecule_info%nmol_types
         WRITE(tab_box_string,'(I4)') 81-10*mc_par%mc_molecule_info%nboxes
         WRITE(tab_string_int,'(I4)') 81-5*mc_par%mc_molecule_info%nmol_types
         format_string="(A,T" // TRIM(ADJUSTL(tab_string)) // "," // TRIM(ADJUSTL(molecule_string)) // "(2X,F8.4))"
         format_box_string="(A,T" // TRIM(ADJUSTL(tab_box_string)) // "," // TRIM(ADJUSTL(box_string)) // "(2X,F8.4))"
         format_string_int="(A,T" // TRIM(ADJUSTL(tab_string)) // "," // TRIM(ADJUSTL(molecule_string)) // "(I3,2X))"
         WRITE ( iw, '( A )' ) ' MC| Monte Carlo Protocol '
         WRITE ( iw, '( A,T71,I10 )' ) ' MC| total number of steps ', &
            mc_par % nstep
         WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| pmvolume ', &
            mc_par % pmvolume
         WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| pmvol_box ', &
            mc_par % pmvol_box
         WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| pmhmc ', &
            mc_par % pmhmc
         WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| pmhmc_box ', &
            mc_par % pmhmc_box
         WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| pmswap ', &
            mc_par % pmswap
         WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| pmavbmc ', &
            mc_par % pmavbmc
         WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| pmtraion ', &
            mc_par % pmtraion
         WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| pmtrans ', &
            mc_par % pmtrans
         WRITE ( iw, '( A,T71,I10 )' ) ' MC| iupvolume ', &
            mc_par % iupvolume
         WRITE ( iw, '( A,T71,I10 )' ) ' MC| nvirial ', &
            mc_par % nvirial
         WRITE ( iw, '( A,T71,I10 )' ) ' MC| iuptrans ', &
            mc_par % iuptrans
         WRITE ( iw, '( A,T71,I10 )' ) ' MC| iprint ', &
            mc_par % iprint
         WRITE ( iw, '( A,T61,A20 )' ) ' MC| ensemble ', &
            ADJUSTR(mc_par % ensemble)
! format string may not be valid if all the molecules are atomic,
! like in HMC
       IF(.NOT. mc_par % lhmc)THEN
            WRITE ( iw, format_string ) ' MC| pmswap_mol ', &
                 mc_par % pmswap_mol
            WRITE ( iw, format_string ) ' MC| pmavbmc_mol ', &
                 mc_par % pmavbmc_mol
            WRITE ( iw, format_string ) ' MC| pbias ', &
                 mc_par % pbias
            WRITE ( iw, format_string ) ' MC| pmtraion_mol ', &
                 mc_par % pmtraion_mol
            WRITE ( iw, format_string ) ' MC| pmtrans_mol ', &
                 mc_par % pmtrans_mol
            WRITE ( iw, format_string ) ' MC| pmrot_mol ', &
                 mc_par % pmrot_mol
            WRITE ( iw, format_string ) ' MC| eta [K]', &
                 mc_par % eta
            WRITE ( iw, format_string ) ' MC| rmbond [angstroms]', &
                 mc_par % rmbond
            WRITE ( iw, format_string ) ' MC| rmangle [degrees]', &
                 mc_par % rmangle
            WRITE ( iw, format_string ) ' MC| rmdihedral [degrees]', &
                 mc_par % rmdihedral
            WRITE ( iw, format_string ) ' MC| rmtrans [angstroms]', &
                 mc_par % rmtrans
            WRITE ( iw, format_string ) ' MC| rmrot [degrees]', &
                 mc_par % rmrot
            WRITE ( iw, format_string_int ) ' MC| AVBMC target atom ', &
                 mc_par % avbmc_atom
            WRITE ( iw, format_string ) ' MC| AVBMC inner cutoff [ang]', &
                 mc_par % avbmc_rmin
            WRITE ( iw, format_string ) ' MC| AVBMC outer cutoff [ang]', &
                 mc_par % avbmc_rmax
         ENDIF
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
         WRITE ( iw, '( A,T71,F10.2 )' ) ' MC| Temperature [K] ', &
            mc_par % temperature
         WRITE ( iw, '( A,T71,F10.5 )' ) ' MC| Pressure [bar] ', &
            mc_par % pressure
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

         WRITE( iw, '(A,I4,A,A)') ' MC| ',mc_par%nswapmoves,&
              ' swap insertions will be attempted ',&
              'per molecular swap move'

         IF(mc_par % lhmc)THEN
            WRITE(iw,'(A)') '********************************************************'
            WRITE(iw,'(A)') '************** ONLY DOING HYBRID MONTE CARLO MOVES **************************'
            WRITE(iw,'(A)') '********************************************************'

         ENDIF
      END IF

! figure out what beta (1/kT) is in atomic units (1/Hartree)
      mc_par % BETA = 1 / mc_par%temperature / boltzmann * joule

      ! 0.9_dp is to avoid overflow or underflow
      mc_par % exp_max_val = 0.9_dp*LOG(HUGE(0.0_dp))
      mc_par % exp_min_val = 0.9_dp*LOG(TINY(0.0_dp))
      mc_par % max_val     = EXP(mc_par % exp_max_val)
      mc_par % min_val     = 0.0_dp

! convert from bar to a.u.
      mc_par%pressure = cp_unit_to_cp2k(mc_par%pressure,"bar",error=error)
! convert from angstrom to a.u.
      DO itype=1,mc_par%mc_molecule_info%nmol_types
! convert from Kelvin to a.u.
!         mc_par % eta = mc_par%eta(itype) * boltzmann / joule
! convert from degrees to radians
         mc_par%rmrot(itype) = mc_par%rmrot(itype)/180.0e0_dp*pi
         mc_par%rmangle(itype) = mc_par%rmangle(itype)/180.0e0_dp*pi
         mc_par%rmdihedral(itype) = mc_par%rmdihedral(itype)/180.0e0_dp*pi
         mc_par%rmtrans(itype) = cp_unit_to_cp2k(mc_par%rmtrans(itype),"angstrom",error=error)
         mc_par%rmbond(itype) = cp_unit_to_cp2k(mc_par%rmbond(itype),"angstrom",error=error)
         mc_par%eta(itype) = cp_unit_to_cp2k(mc_par%eta(itype),"K_e",error=error)
         mc_par%avbmc_rmin(itype) = cp_unit_to_cp2k(mc_par%avbmc_rmin(itype),"angstrom",error=error)
         mc_par%avbmc_rmax(itype) = cp_unit_to_cp2k(mc_par%avbmc_rmax(itype),"angstrom",error=error)
      ENDDO
      mc_par%rmvolume = cp_unit_to_cp2k(mc_par%rmvolume,"angstrom^3",error=error)
      mc_par%discrete_step = cp_unit_to_cp2k(mc_par%discrete_step,"angstrom",error=error)

! end the timing
      CALL timestop(handle)

END SUBROUTINE read_mc_section

! *****************************************************************************
!> \brief finds the largest interaction cutoff value in a classical simulation
!>      so we know the smallest size we can make the box in a volume move
!> \param mc_par the structure that will store the parameters
!> \param force_env the force environment that we'll grab the rcut parameter
!>                   out of
!> \param lterminate set to .TRUE. if one of the sides of the box is
!>            less than twice the cutoff
!>
!>    Suitable for parallel.
!> \author MJM
! *****************************************************************************
SUBROUTINE find_mc_rcut ( mc_par, force_env, lterminate )

    TYPE(mc_simpar_type), INTENT(INOUT)      :: mc_par
    TYPE(force_env_type), POINTER            :: force_env
    LOGICAL, INTENT(OUT)                     :: lterminate

    CHARACTER(len=*), PARAMETER :: routineN = 'find_mc_rcut', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: itype, jtype
    REAL(KIND=dp)                            :: rcutsq_max
    REAL(KIND=dp), DIMENSION(1:3)            :: abc
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type)                      :: error
    TYPE(fist_environment_type), POINTER     :: fist_env
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(pair_potential_pp_type), POINTER    :: potparm

      NULLIFY(cell,potparm,fist_nonbond_env,fist_env)

      lterminate=.FALSE.
      CALL force_env_get(force_env,fist_env=fist_env,error=error)
      CALL fist_env_get(fist_env,cell=cell,&
         fist_nonbond_env=fist_nonbond_env,error=error)
      CALL fist_nonbond_env_get (fist_nonbond_env, potparm=potparm,error=error)
      CALL get_cell(cell,abc=abc)

! find the largest value of rcutsq
      rcutsq_max=0.0e0_dp
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

! *****************************************************************************
!> \brief figures out the number of total molecules, the number of atoms in each
!>      molecule, an array with the molecule types, etc...a lot of information
!>      that we need.  I did this because I use multiple force_env (simulation
!>      boxes) for MC, and they don't know about each other.
!> \param force_env the pointer containing all the force environments in the
!>        simulation
!> \param mc_molecule_info the structure that will hold all the information
!>          for the molecule types
!>
!>    Suitable for parallel.
!> \author MJM
! *****************************************************************************
SUBROUTINE mc_determine_molecule_info(force_env,mc_molecule_info,error,box_number,&
     coordinates_empty)

    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env
    TYPE(mc_molecule_info_type), POINTER     :: mc_molecule_info
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, INTENT(IN), OPTIONAL            :: box_number
    REAL(dp), DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: coordinates_empty

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mc_determine_molecule_info', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name
    CHARACTER(LEN=default_string_length), &
      ALLOCATABLE, DIMENSION(:)              :: names_init
    INTEGER :: first_mol, iatom, ibox, imol, imolecule, istat, itype, &
      iunique, iunit, jtype, last_mol, natom, natoms_large, nbend, nbond, &
      nboxes, nmol_types, nmolecule, ntorsion, ntypes, skip_box, &
      this_molecule, total
    INTEGER, DIMENSION(:), POINTER           :: molecule_list
    LOGICAL                                  :: lnew_type
    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(mol_kind_new_list_p_type), &
      DIMENSION(:), POINTER                  :: molecule_kinds_new
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(particle_list_p_type), &
      DIMENSION(:), POINTER                  :: particles

! how many simulation boxes do we have?

    nboxes=SIZE(force_env(:))

! allocate the pointer
    ALLOCATE(mc_molecule_info)
    mc_molecule_info%nboxes=nboxes

! if box_number is present, that box is supposed to be empty,
! so we don't count it in anything
    IF(PRESENT(box_number)) THEN
       skip_box=box_number
    ELSE
       skip_box=0
    ENDIF

! we need to figure out how many different kinds of molecules we have...
! the fun stuff comes in when box 1 is missing a molecule type...I'll
! distinguish molecules based on names from the .psf files...
! this is only really a problem when we have more than one simulation
! box
    NULLIFY(subsys,molecule_kinds_new,molecule_kind)

    ALLOCATE(particles(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "particles",nboxes)
    ALLOCATE (molecule_kinds_new(1:nboxes),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "molecule_kinds_new",nboxes)

 ! find out how many types we have
    ntypes=0
    DO ibox=1,nboxes
       IF(ibox == skip_box) CYCLE
       CALL force_env_get(force_env(ibox)%force_env,&
            subsys=subsys,error=error)
       CALL cp_subsys_get(subsys, &
            molecule_kinds_new=molecule_kinds_new(ibox)%list,&
            particles=particles(ibox)%list,error=error)
       ntypes=ntypes+SIZE(molecule_kinds_new(ibox)%list%els(:))
    ENDDO

    ALLOCATE (names_init(1:ntypes),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "names_init",ntypes*default_string_length)

! now get the names of all types so we can figure out how many distinct
! types we have
    itype=1
    natoms_large=0
    DO ibox=1,nboxes
       IF(ibox == skip_box) CYCLE
       DO imolecule=1,SIZE(molecule_kinds_new(ibox)%list%els(:))
          molecule_kind=> molecule_kinds_new(ibox)%list%els(imolecule)
          CALL get_molecule_kind(molecule_kind,name=names_init(itype),&
               natom=natom)
          IF(natom .GT. natoms_large) natoms_large=natom
          itype=itype+1
       ENDDO
    ENDDO

    nmol_types=0
    DO itype=1,ntypes
       lnew_type=.TRUE.
       DO jtype=1,itype-1
          IF(TRIM(names_init(itype)) .EQ. TRIM(names_init(jtype)))&
               lnew_type=.FALSE.
       ENDDO
       IF(lnew_type) THEN
          nmol_types=nmol_types+1
       ELSE
          names_init(itype)=''
       ENDIF
    ENDDO

! allocate arrays for the names of the molecules, the number of atoms, and
! the conformational probabilities
    mc_molecule_info%nmol_types=nmol_types
    ALLOCATE(mc_molecule_info%names(1:nmol_types),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_molecule_info%names",nmol_types*default_string_length)
    ALLOCATE(mc_molecule_info%nunits(1:nmol_types),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_molecule_info%nunits",nmol_types*int_size)
    ALLOCATE(mc_molecule_info%conf_prob(1:3,1:nmol_types),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_molecule_info%conf_prob",3*nmol_types*dp_size)
    ALLOCATE(mc_molecule_info%nchains(1:nmol_types,1:nboxes),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_molecule_info%nchains",nboxes*nmol_types*dp_size)
    ALLOCATE(mc_molecule_info%nunits_tot(1:nboxes),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_molecule_info%nunits_tot",nboxes*int_size)
    ALLOCATE(mc_molecule_info%atom_names(1:natoms_large,1:nmol_types),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_molecule_info%atom_names",natoms_large*nmol_types*default_string_length)
    ALLOCATE(mc_molecule_info%mass(1:natoms_large,1:nmol_types),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_molecule_info%mass",natoms_large*nmol_types*dp_size)

! now record all the information for every molecule type
    iunique=0
    itype=0
    DO ibox=1,nboxes
       IF(ibox == skip_box) CYCLE
       DO imolecule=1,SIZE(molecule_kinds_new(ibox)%list%els(:))
          itype=itype+1
          IF(names_init(itype) .EQ. '') CYCLE
          iunique=iunique+1
          mc_molecule_info%names(iunique)=names_init(itype)
          molecule_kind=> molecule_kinds_new(ibox)%list%els(imolecule)
          CALL get_molecule_kind(molecule_kind,natom=mc_molecule_info%nunits(iunique),&
               nbond=nbond,nbend=nbend,ntorsion=ntorsion,atom_list=atom_list)

          DO iatom=1,mc_molecule_info%nunits(iunique)
             atomic_kind => atom_list(iatom)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  name=mc_molecule_info%atom_names(iatom,iunique),&
                  mass=mc_molecule_info%mass(iatom,iunique))
          ENDDO

! compute the probabilities of doing any particular kind of conformation change
          total=nbond+nbend+ntorsion

          IF(total == 0) THEN
             mc_molecule_info%conf_prob(:,iunique)=0.0e0_dp
          ELSE
             mc_molecule_info%conf_prob(1,iunique)=REAL(nbond,dp)/REAL(total,dp)
             mc_molecule_info%conf_prob(2,iunique)=REAL(nbend,dp)/REAL(total,dp)
             mc_molecule_info%conf_prob(3,iunique)=REAL(ntorsion,dp)/REAL(total,dp)
          ENDIF

       ENDDO
    ENDDO

! figure out the empty coordinates
    IF(PRESENT(coordinates_empty)) THEN
       ALLOCATE(coordinates_empty(1:3,1:mc_molecule_info%nunits(1)),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "coordinates_empty",3*mc_molecule_info%nunits(1)*dp_size)
       DO iunit=1,mc_molecule_info%nunits(1)
          coordinates_empty(1:3,iunit)=particles(1)%list%els(iunit)%r(1:3)
       ENDDO
    ENDIF

! now we need to figure out the total number of molecules of each kind in
! each box, and the total number of interaction sites in each box
    mc_molecule_info%nchains(:,:)=0
    DO iunique=1,nmol_types
       DO ibox=1,nboxes
          IF(ibox == skip_box) CYCLE
          DO imolecule=1,SIZE(molecule_kinds_new(ibox)%list%els(:))
             molecule_kind=> molecule_kinds_new(ibox)%list%els(imolecule)
             CALL get_molecule_kind(molecule_kind,nmolecule=nmolecule,&
                  name=name)
             IF(TRIM(name) .NE. mc_molecule_info%names(iunique) ) CYCLE
             mc_molecule_info%nchains(iunique,ibox)=mc_molecule_info%nchains(iunique,ibox)+nmolecule
          ENDDO
       ENDDO
    ENDDO
    mc_molecule_info%nchain_total=0
    DO ibox=1,nboxes
       mc_molecule_info%nunits_tot(ibox)=0
       IF(ibox == skip_box) CYCLE
       DO iunique=1,nmol_types
          mc_molecule_info%nunits_tot(ibox)=mc_molecule_info%nunits_tot(ibox)+&
               mc_molecule_info%nunits(iunique)*mc_molecule_info%nchains(iunique,ibox)
       ENDDO
       mc_molecule_info%nchain_total=mc_molecule_info%nchain_total+SUM(mc_molecule_info%nchains(:,ibox))
    ENDDO

! now we need to figure out which type every molecule is,
! and which force_env (box) it's in
    ALLOCATE(mc_molecule_info%mol_type(1:mc_molecule_info%nchain_total),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mol_type",mc_molecule_info%nchain_total*int_size)
    ALLOCATE(mc_molecule_info%in_box(1:mc_molecule_info%nchain_total),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "in_box",mc_molecule_info%nchain_total*int_size)

    last_mol=0
    DO ibox=1,nboxes
       IF(ibox == skip_box) CYCLE
       first_mol=last_mol+1
       last_mol=first_mol+SUM(mc_molecule_info%nchains(:,ibox))-1
       mc_molecule_info%in_box(first_mol:last_mol)=ibox
       DO imolecule=1,SIZE(molecule_kinds_new(ibox)%list%els(:))
          molecule_kind=> molecule_kinds_new(ibox)%list%els(imolecule)
          CALL get_molecule_kind(molecule_kind,nmolecule=nmolecule,&
               name=name,molecule_list=molecule_list)
! figure out which molecule number this is
          this_molecule=0
          DO iunique=1,nmol_types
             IF(TRIM(name) .EQ. TRIM(mc_molecule_info%names(iunique))) THEN
                this_molecule=iunique
             ENDIF
          ENDDO
          DO imol=1,SIZE(molecule_list(:))
             mc_molecule_info%mol_type(first_mol+molecule_list(imol)-1)=this_molecule
          ENDDO
       ENDDO
    ENDDO

! get rid of stuff we don't need
    DEALLOCATE(names_init,stat=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"names_init")
    DEALLOCATE(molecule_kinds_new,stat=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"molecule_kinds_new")
    DEALLOCATE(particles,stat=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"particles")

END SUBROUTINE MC_DETERMINE_MOLECULE_INFO

! *****************************************************************************
!> \brief deallocates all the arrays in the mc_molecule_info_type
!> \param mc_molecule_info the structure we wish to deallocate
!>
!>    Suitable for parallel.
!> \author MJM
! *****************************************************************************
SUBROUTINE mc_molecule_info_destroy(mc_molecule_info)

    TYPE(mc_molecule_info_type), POINTER     :: mc_molecule_info

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mc_molecule_info_destroy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    DEALLOCATE(mc_molecule_info%nchains,stat=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"nchains")
    DEALLOCATE(mc_molecule_info%nunits,stat=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"nunits")
    DEALLOCATE(mc_molecule_info%mol_type,stat=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"mol_type")
    DEALLOCATE(mc_molecule_info%in_box,stat=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"in_box")
    DEALLOCATE(mc_molecule_info%names,stat=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"names")
    DEALLOCATE(mc_molecule_info%atom_names,stat=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"atom_names")
    DEALLOCATE(mc_molecule_info%conf_prob,stat=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"conf_prob")
    DEALLOCATE(mc_molecule_info%nunits_tot,stat=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"nunits_tot")
    DEALLOCATE(mc_molecule_info%mass,stat=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"mass")

    DEALLOCATE(mc_molecule_info,stat=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"mc_molecule_info")
    NULLIFY(mc_molecule_info)

  END SUBROUTINE mc_molecule_info_destroy

! *****************************************************************************
  SUBROUTINE mc_input_parameters_check ( mc_par)

! *****************************************************************************
!> \brief accesses the private elements of the mc_molecule_info_type
!> \param mc_molecule_info the structure you want the parameters for
!> \param nmol_types the number of molecule types in the simulation
!>
!>    Suitable for parallel.
!> \author MJM
    TYPE(mc_simpar_type), POINTER            :: mc_par

    CHARACTER(len=*), PARAMETER :: routineN = 'mc_input_parameters_check', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imol_type, nboxes, &
                                                nchain_total, nmol_types
    INTEGER, DIMENSION(:), POINTER           :: nunits
    LOGICAL                                  :: lhmc
    TYPE(mc_molecule_info_type), POINTER     :: mc_molecule_info

    CALL get_mc_par(mc_par,mc_molecule_info=mc_molecule_info,lhmc=lhmc)

    CALL get_mc_molecule_info(mc_molecule_info,nmol_types=nmol_types,&
         nboxes=nboxes,nunits=nunits,nchain_total=nchain_total)

! if we're only doing HMC, we can skip these checks
    IF(lhmc) RETURN

! the final value of these arrays needs to be 1.0, otherwise there is
! a chance we won't select a molecule type for a move
    IF(mc_par%pmavbmc_mol(nmol_types) .LT. 0.9999E0_dp) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
            'The last value of PMAVBMC_MOL needs to be 1.0')
    ENDIF
    IF(mc_par%pmswap_mol(nmol_types) .LT. 0.9999E0_dp) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
            'The last value of PMSWAP_MOL needs to be 1.0')
    ENDIF
    IF(mc_par%pmtraion_mol(nmol_types) .LT. 0.9999E0_dp) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
            'The last value of PMTRAION_MOL needs to be 1.0')
    ENDIF
    IF(mc_par%pmtrans_mol(nmol_types) .LT. 0.9999E0_dp) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
            'The last value of PMTRANS_MOL needs to be 1.0')
    ENDIF
    IF(mc_par%pmrot_mol(nmol_types) .LT. 0.9999E0_dp) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
            'The last value of PMROT_MOL needs to be 1.0')
    ENDIF

! check to make sure we have all the desired values for various

 ! some ensembles require multiple boxes or molecule types
    SELECT CASE(mc_par%ensemble)
    CASE("GEMC_NPT")
       IF(nmol_types .LE. 1) &
            CALL stop_program(routineN,moduleN,__LINE__,&
            'Cannot have GEMC-NPT simulation with only one molecule type')
       IF(nboxes .LE. 1) &
            CALL stop_program(routineN,moduleN,__LINE__,&
            'Cannot have GEMC-NPT simulation with only one box')
    CASE("GEMC_NVT")
       IF(nboxes .LE. 1) &
            CALL stop_program(routineN,moduleN,__LINE__,&
            'Cannot have GEMC-NVT simulation with only one box')
    CASE("TRADITIONAL")
       IF(mc_par%pmswap .GT. 0.0E0_dp)&
            CALL stop_program(routineN,moduleN,__LINE__,&
            'You cannot do swap moves in a system with only one box')
    CASE("VIRIAL")
       IF(nchain_total .NE. 2)&
            CALL stop_program(routineN,moduleN,__LINE__,&
            'You need exactly two molecules in the box to compute the second virial.')
    END SELECT

! can't choose an AVBMC target atom number higher than the number
! of units in that molecule
    DO imol_type=1,nmol_types
       IF(mc_par%avbmc_atom(imol_type) .GT. nunits(imol_type)) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
               'Cannot have avbmc_atom higher than the number of atoms for that molecule type!')
       ENDIF
    ENDDO

  END SUBROUTINE mc_input_parameters_check

END MODULE mc_types

