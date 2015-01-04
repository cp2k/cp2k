!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      10.2005 split input_cp2k into smaller modules [fawzi]
!> \author teo & fawzi
! *****************************************************************************
MODULE input_cp2k_motion
  USE bibliography,                    ONLY: Byrd1995,&
                                             Henkelman1999
  USE cp_output_handling,              ONLY: add_last_numeric,&
                                             cp_print_key_section_create,&
                                             debug_print_level,&
                                             high_print_level,&
                                             low_print_level,&
                                             medium_print_level
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE input_constants,                 ONLY: &
       default_bfgs_method_id, default_cell_direct_id, &
       default_cell_geo_opt_id, default_cell_md_id, default_cg_method_id, &
       default_dimer_method_id, default_lbfgs_method_id, &
       default_minimization_method_id, default_ts_method_id, do_mc_gemc_npt, &
       do_mc_gemc_nvt, do_mc_traditional, do_mc_virial, fmt_id_pdb, &
       fmt_id_xyz, gaussian, helium_cell_shape_cube, &
       helium_cell_shape_octahedron, ls_2pnt, ls_3pnt, ls_fit, ls_gold, &
       ls_none, numerical, transformation_normal, transformation_stage
  USE input_cp2k_constraints,          ONLY: create_constraint_section
  USE input_cp2k_free_energy,          ONLY: create_fe_section
  USE input_cp2k_md,                   ONLY: create_md_section
  USE input_cp2k_motion_print,         ONLY: add_format_keyword,&
                                             create_motion_print_section
  USE input_cp2k_neb,                  ONLY: create_band_section
  USE input_cp2k_thermostats,          ONLY: create_coord_section,&
                                             create_gle_section,&
                                             create_velocity_section
  USE input_cp2k_tmc,                  ONLY: create_TMC_section
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: integer_t,&
                                             real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_motion'

PUBLIC :: create_motion_section

CONTAINS

! *****************************************************************************
!> \brief creates the motion section
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_motion_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_motion_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="motion",&
            description="This section defines a set of tool connected with the motion of the nuclei.",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection)

       CALL create_geoopt_section(subsection,label="GEO_OPT",&
            description="This section sets the environment of the geometry optimizer.",&
            just_optimizers=.FALSE.,error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_cell_opt_section(subsection, error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_shellcore_opt_section(subsection, error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_md_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_driver_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_fe_section(subsection,error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_constraint_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_fp_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_mc_section(subsection, error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_TMC_section(subsection, error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_pint_section(subsection, error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_band_section(subsection, error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_motion_print_section(subsection, error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
    END IF

  END SUBROUTINE create_motion_section

! *****************************************************************************
!> \brief creates the Monte Carlo section
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author matt
! *****************************************************************************
  SUBROUTINE create_mc_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mc_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="mc",&
            description="This section sets parameters to set up a MonteCarlo calculation.",&
            n_keywords=10, n_subsections=2, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword,subsection)

       CALL keyword_create(keyword, name="NSTEP",&
            description="Specifies the number of MC cycles.",&
            usage="NSTEP {integer}",&
            default_i_val=100,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="IPRINT",&
            description="Prints coordinate/cell/etc information every IPRINT steps.",&
            usage="IPRINT {integer}",&
            default_i_val=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NMOVES",&
            description="Specifies the number of classical moves between energy evaluations. ",&
            usage="NMOVES {integer}",&
            default_i_val=4,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NSWAPMOVES",&
            description="How many insertions to try per swap move.",&
            usage="NSWAPMOVES {integer}",&
            default_i_val=16,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LBIAS",&
            description="Dictates if we presample moves with a different potential.",&
            usage="LBIAS {logical}",&
            default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LSTOP",&
            description="Makes nstep in terms of steps, instead of cycles.",&
            usage="LSTOP {logical}",&
            default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LDISCRETE",&
            description="Changes the volume of the box in discrete steps, one side at a time.",&
            usage="LDISCRETE {logical}",&
            default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART",&
            description="Read initial configuration from restart file.",&
            usage="RESTART {logical}",&
            default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NVIRIAL",&
            description="Use this many random orientations to compute the second virial coefficient (ENSEMBLE=VIRIAL)",&
            usage="NVIRIAL {integer}",&
            default_i_val=1000,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ENSEMBLE",&
            description="Specify the type of simulation",&
            usage="PROGRAM (TRADITIONAL|GEMC_NVT|GEMC_NPT|VIRIAL)",&
            enum_c_vals=s2a( "TRADITIONAL","GEMC_NVT","GEMC_NPT","VIRIAL"),&
            enum_i_vals=(/do_mc_traditional,do_mc_gemc_nvt,do_mc_gemc_npt,do_mc_virial/),&
            default_i_val=do_mc_traditional,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART_FILE_NAME",&
            description="Name of the restart file for MC information.",&
            usage="RESTART_FILE_NAME {filename}",&
            default_lc_val="",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MOVES_FILE_NAME",&
            description="The file to print the move statistics to.",&
            usage="MOVES_FILE_NAME {filename}",&
            default_lc_val="",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MOLECULES_FILE_NAME",&
            description="The file to print the number of molecules to.",&
            usage="MOLECULES_FILE_NAME {filename}",&
            default_lc_val="",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="COORDINATE_FILE_NAME",&
            description="The file to print the current coordinates to.",&
            usage="COORDINATE_FILE_NAME {filename}",&
            default_lc_val="",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ENERGY_FILE_NAME",&
            description="The file to print current energies to.",&
            usage="ENERGY_FILE_NAME {filename}",&
            default_lc_val="",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DATA_FILE_NAME",&
            description="The file to print current configurational info to.",&
            usage="DATA_FILE_NAME {filename}",&
            default_lc_val="",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CELL_FILE_NAME",&
            description="The file to print current cell length info to.",&
            usage="CELL_FILE_NAME {filename}",&
            default_lc_val="",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_DISP_FILE_NAME",&
            description="The file to print current maximum displacement info to.",&
            usage="MAX_DISP_FILE_NAME {filename}",&
            default_lc_val="",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BOX2_FILE_NAME",&
            description="For GEMC, the name of the input file for the other box.",&
            usage="BOX2_FILE_NAME {filename}",&
            default_lc_val="",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PRESSURE",&
            description="The pressure for NpT simulations, in bar.",&
            usage="PRESSURE {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TEMPERATURE",&
            description="The temperature of the simulation, in Kelvin.",&
            usage="TEMPERATURE {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="VIRIAL_TEMPS",&
            description="The temperatures you wish to compute the virial coefficient for.  Only used if ensemble=VIRIAL.",&
            usage="VIRIAL_TEMPS {real} {real} ... ",&
            n_var=-1,type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DISCRETE_STEP",&
            description="The size of the discrete volume move step, in angstroms.",&
            usage="DISCRETE_STEP {real}",&
            default_r_val=1.0E0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ETA",&
            description="The free energy bias (in Kelvin) for swapping a molecule of each type into this box.",&
            usage="ETA {real} {real} ... ",&
            n_var=-1,type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RANDOMTOSKIP",&
            description="Number of random numbers from the acceptance/rejection stream to skip",&
            usage="RANDOMTOSKIP {integer}",&
            default_i_val=0,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_avbmc_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_move_prob_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_update_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_max_disp_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

     END IF

  END SUBROUTINE create_mc_section

! *****************************************************************************
!> \brief ...
!> \param section will contain the AVBMC parameters for MC
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author matt
! *****************************************************************************
  SUBROUTINE create_avbmc_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_avbmc_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="avbmc",&
            description="Parameters for Aggregation Volume Bias Monte Carlo (AVBMC) "//&
            "which explores cluster formation and destruction. "//&
            "Chen and Siepmann, J. Phys. Chem. B 105, 11275-11282 (2001).",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword )

       CALL keyword_create(keyword, name="PBIAS",&
            description="The probability of swapping to an inner region in an AVBMC swap move for each molecule type.",&
            usage="PBIAS {real} {real} ... ",&
            n_var=-1,type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVBMC_ATOM",&
            description="The target atom for an AVBMC swap move for each molecule type.",&
            usage="AVBMC_ATOM {integer} {integer} ... ",&
            n_var=-1,type_of_var=integer_t,required=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVBMC_RMIN",&
            description="The inner radius for an AVBMC swap move, in angstroms for every molecule type.",&
            usage="AVBMC_RMIN {real} {real} ... ",&
            n_var=-1,type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVBMC_RMAX",&
            description="The outer radius for an AVBMC swap move, in angstroms, for every molecule type.",&
            usage="AVBMC_RMAX {real} {real} ... ",&
            n_var=-1,type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_avbmc_section

! *****************************************************************************
!> \brief ...
!> \param section will contain the probabilities for attempting each move
!>        type in Monte Carlo
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author matt
! *****************************************************************************
  SUBROUTINE create_move_prob_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_move_prob_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="move_probabilities",&
            description="Parameters for fraction of moves performed for each move type.",&
            n_keywords=5, n_subsections=2, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword,subsection )

       CALL keyword_create(keyword, name="PMHMC",&
            description="The probability of attempting a hybrid MC move.",&
            usage="PMHMC {real}",&
            type_of_var=real_t, default_r_val=0.0E0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PMTRANS",&
            description="The probability of attempting a molecule translation.",&
            usage="PMTRANS {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PMAVBMC",&
            description="The probability of attempting an AVBMC swap move.",&
            usage="PMAVBMC {real}",&
            default_r_val=0.0E0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PMTRAION",&
            description="The probability of attempting a conformational change.",&
            usage="PMTRAION {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PMSWAP",&
            description="The probability of attempting a swap move.",&
            usage="PMSWAP {real}",&
            type_of_var=real_t, default_r_val=0.0E0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PMVOLUME",&
            description="The probability of attempting a volume move.",&
            usage="PMVOLUME {real}",&
            type_of_var=real_t, default_r_val=0.0E0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_mol_prob_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_box_prob_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_move_prob_section

! *****************************************************************************
!> \brief ...
!> \param section will contain the probabilities for attempting various moves
!>        on the various molecule types present in the system
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author matt
! *****************************************************************************
  SUBROUTINE create_mol_prob_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mol_prob_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="mol_probabilities",&
            description="Probabilities of attempting various moves types on "//&
            "the various molecular types present in the simulation.",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword )

       CALL keyword_create(keyword, name="PMAVBMC_MOL",&
            description="The probability of attempting an AVBMC swap move on each molecule type.",&
            usage="PMAVBMC_MOL {real} {real} ... ",&
            n_var=-1,type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PMSWAP_MOL",&
            description="The probability of attempting a molecule swap of a given molecule type.",&
            usage="PMSWAP_MOL {real} {real} ... ",&
            n_var=-1,type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PMROT_MOL",&
            description="The probability of attempting a molecule rotation of a given molecule type.",&
            usage="PMROT_MOL {real} {real} ... ",&
            n_var=-1,type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PMTRAION_MOL",&
            description="The probability of attempting a conformational change of a given molecule type.",&
            usage="PMTRAION_MOL {real} {real} ... ",&
            n_var=-1,type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PMTRANS_MOL",&
            description="The probability of attempting a molecule translation of a given molecule type.",&
            usage="PMTRANS_MOL {real} {real} ... ",&
            n_var=-1,type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_mol_prob_section

! *****************************************************************************
!> \brief ...
!> \param section will contain the probabilities for attempting various moves
!>        on the box where the variable is present
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author matt
! *****************************************************************************
  SUBROUTINE create_box_prob_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_box_prob_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="BOX_PROBABILITIES",&
            description="Probabilities of attempting various moves types on "//&
            "the box.",&
            n_keywords=2, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword )

       CALL keyword_create(keyword, name="PMHMC_BOX",&
            description="The probability of attempting a HMC move on this box.",&
            usage="PMHMC_BOX {real}",&
            type_of_var=real_t, default_r_val=1.0E0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PMVOL_BOX",&
            description="The probability of attempting a volume move on this box (GEMC_NpT).",&
            usage="PMVOL_BOX {real}",&
            type_of_var=real_t, default_r_val=1.0E0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


    END IF
  END SUBROUTINE create_box_prob_section

! *****************************************************************************
!> \brief ...
!> \param section will contain the frequency for updating maximum
!>        displacements for various moves
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author matt
! *****************************************************************************
  SUBROUTINE create_update_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_update_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="MOVE_UPDATES",&
            description="Frequency for updating move maximum displacements.",&
            n_keywords=2, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword )

       CALL keyword_create(keyword, name="IUPVOLUME",&
            description="Every iupvolume steps update maximum volume displacement.",&
            usage="IUPVOLUME {integer}",&
            default_i_val=10000,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="IUPTRANS",&
            description="Every iuptrans steps update maximum translation/rotation/" //&
            "configurational changes.",&
            usage="IUPTRANS {integer}",&
            default_i_val=10000,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


    END IF
  END SUBROUTINE create_update_section

! *****************************************************************************
!> \brief ...
!> \param section will contain the maximum displacements for various moves
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author matt
! *****************************************************************************
  SUBROUTINE create_max_disp_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_max_disp_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="max_displacements",&
            description="The maximum displacements for all attempted moves.",&
            n_keywords=1, n_subsections=2, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection)

       CALL create_mol_disp_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_box_disp_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_max_disp_section

! *****************************************************************************
!> \brief ...
!> \param section will contain the maximum displacements for all moves which
!>        require a value for each molecule type
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author matt
! *****************************************************************************
  SUBROUTINE create_mol_disp_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mol_disp_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="mol_displacements",&
            description="Maximum displacements for every move type that requires "//&
            "a value for each molecular type in the simulation.",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword )

       CALL keyword_create(keyword, name="RMBOND",&
            description="Maximum bond length displacement, in angstroms, for each molecule type.",&
            usage="RMBOND {real} {real} ... ",&
            n_var=-1,type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMANGLE",&
            description="Maximum bond angle displacement, in degrees, for each molecule type.",&
            usage="RMANGLE {real} {real} ...",&
            n_var=-1,type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMDIHEDRAL",&
            description="Maximum dihedral angle distplacement, in degrees, for each molecule type.",&
            usage="RMDIHEDRAL {real} {real} ... ",&
            n_var=-1,type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMROT",&
            description="Maximum rotational displacement, in degrees, for each molecule type.",&
            usage="RMROT {real} {real} ... ",&
            n_var=-1,type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMTRANS",&
            description="Maximum translational displacement, in angstroms, for each molecule type.",&
            usage="RMTRANS {real} {real} ...",&
            n_var=-1,type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_mol_disp_section

! *****************************************************************************
!> \brief ...
!> \param section will contain the maximum displacements for any move that
!>        is done on each simulation box
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author matt
! *****************************************************************************
  SUBROUTINE create_box_disp_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_box_disp_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="BOX_DISPLACEMENTS",&
            description="Maximum displacements for any move that is performed on each" // &
            " simulation box.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword )

       CALL keyword_create(keyword, name="RMVOLUME",&
            description="Maximum volume displacement, in angstrom**3.",&
            usage="RMVOLUME {real}",&
            type_of_var=real_t, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


    END IF
  END SUBROUTINE create_box_disp_section

! *****************************************************************************
!> \brief creates the geometry optimization section
!> \param section the section to be created
!> \param label ...
!> \param description ...
!> \param just_optimizers ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  RECURSIVE SUBROUTINE create_geoopt_section(section,label,description,just_optimizers,error)
    TYPE(section_type), POINTER              :: section
    CHARACTER(LEN=*), INTENT(IN)             :: label, description
    LOGICAL, INTENT(IN)                      :: just_optimizers
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_geoopt_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name=label, description=description,&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)
       IF (.NOT.just_optimizers) THEN
          CALL keyword_create(keyword, name="TYPE",&
               description="Specify which kind of geometry optimization to perform",&
               usage="TYPE (MINIMIZATION|TRANSITION_STATE)",&
               enum_c_vals=s2a( "MINIMIZATION","TRANSITION_STATE"),&
               enum_desc=s2a("Performs a geometry minimization.",&
                             "Performs a transition state optimization."),&
               enum_i_vals=(/default_minimization_method_id,default_ts_method_id/),&
               default_i_val=default_minimization_method_id,error=error)
          CALL section_add_keyword(section,keyword,error=error)
          CALL keyword_release(keyword,error=error)
       END IF

       CALL keyword_create(keyword, name="OPTIMIZER",&
            variants=(/"MINIMIZER"/),&
            citations=(/Byrd1995/),&
            description="Specify which method to use to perform a geometry optimization.",&
            usage="OPTIMIZER {BFGS|LBFGS|CG}",&
            enum_c_vals=s2a( "BFGS","LBFGS","CG"),&
            enum_desc=s2a("Most efficient minimizer, but only for 'small' systems, "//&
            "as it relies on diagonalization of a full Hessian matrix",&
            "Limit memory variant of the above, can also be applied to large systems, not as well fine-tuned",&
            "conjugate gradients, robust minimizer (depending on the line search) also OK for large systems"),&
            enum_i_vals=(/default_bfgs_method_id,default_lbfgs_method_id,default_cg_method_id/),&
            default_i_val=default_bfgs_method_id,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_ITER",&
            description="Specifies the maximum number of geometry optimization steps. "//&
            "One step might imply several force evaluations for the CG and LBFGS optimizers.",&
            usage="MAX_ITER {integer}",&
            default_i_val=200,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_DR",&
            description="Convergence criterium for the maximum geometry change "//&
            "between the current and the last optimizer iteration.",&
            usage="MAX_DR {real}",&
            default_r_val=0.0030_dp,unit_str="bohr",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_FORCE",&
            description="Convergence criterium for the maximum force component of the current configuration.",&
            usage="MAX_FORCE {real}",&
            default_r_val=0.00045_dp,unit_str="hartree/bohr",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMS_DR",&
            description="Convergence criterium for the root mean square (RMS) geometry"//&
            " change between the current and the last optimizer iteration.",&
            usage="RMS_DR {real}",unit_str="bohr",&
            default_r_val=0.0015_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMS_FORCE",&
            description="Convergence criterium for the root mean square (RMS) force of the current configuration.",&
            usage="RMS_FORCE {real}",unit_str="hartree/bohr",&
            default_r_val=0.00030_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,name="step_start_val",&
            description="The starting step value for the "//TRIM(label)//" module.",&
            usage="step_start_val <integer>",default_i_val=0,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_lbfgs_section(subsection, error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_cg_section(subsection, error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_bfgs_section(subsection, error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       IF (.NOT.just_optimizers) THEN
          ! Transition states section
          CALL create_ts_section(subsection, error)
          CALL section_add_subsection(section,subsection,error=error)
          CALL section_release(subsection,error=error)

          ! Create the PRINT subsection
          NULLIFY(subsection)
          CALL section_create(subsection,name="PRINT",&
               description="Controls the printing properties during a geometry optimization run",&
               n_keywords=0, n_subsections=1, repeats=.TRUE., required=.FALSE.,&
               error=error)
          NULLIFY(print_key)
          CALL cp_print_key_section_create(print_key,"program_run_info",&
               description="Controls the printing of basic information during the Geometry Optimization", &
               print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
               error=error)
          CALL section_add_subsection(subsection,print_key,error=error)
          CALL section_release(print_key,error=error)
          CALL section_add_subsection(section,subsection,error=error)
          CALL section_release(subsection,error=error)
       END IF

    END IF
  END SUBROUTINE create_geoopt_section

! *****************************************************************************
!> \brief creates the section for the shell-core optimization
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Caino
! *****************************************************************************
  SUBROUTINE create_shellcore_opt_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_shellcore_opt_section', &
      routineP = moduleN//':'//routineN

    TYPE(section_type), POINTER              :: print_key, subsection

    CALL create_geoopt_section(section, label="SHELL_OPT",&
         description="This section sets the environment for the optimization of the shell-core distances"//&
         " that might turn to be necessary along a MD run using a shell-model potential. "//&
         " The optimization procedure is activated when at least one of the shell-core "//&
         "pairs becomes too elongated,  i.e. when the assumption of point dipole is not longer valid.",&
         just_optimizers=.TRUE.,error=error)

    NULLIFY(print_key, subsection)


    ! Create the PRINT subsection
    NULLIFY(subsection)
    CALL section_create(subsection,name="PRINT",&
         description="Controls the printing properties during a shell-core optimization procedure",&
         n_keywords=0, n_subsections=1, repeats=.TRUE., required=.FALSE.,&
         error=error)
    NULLIFY(print_key)
    CALL cp_print_key_section_create(print_key,"program_run_info",&
         description="Controls the printing of basic information during the Optimization", &
         print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
         error=error)
    CALL section_add_subsection(subsection,print_key,error=error)
    CALL section_release(print_key,error=error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

  END SUBROUTINE create_shellcore_opt_section

! *****************************************************************************
!> \brief creates the section for the cell optimization
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] - University of Zurich - 03.2008
! *****************************************************************************
  SUBROUTINE create_cell_opt_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_cell_opt_section', &
      routineP = moduleN//':'//routineN

    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    CALL create_geoopt_section(section, label="CELL_OPT",&
         description="This section sets the environment for the optimization of the simulation cell."//&
         " Two possible schemes are available: (1) Zero temperature optimization; "//&
         " (2) Finite temperature optimization. ",just_optimizers=.TRUE.,error=error)

    NULLIFY(keyword,print_key,subsection)
    CALL keyword_create(keyword, name="TYPE",&
         description="Specify which kind of method to use for the optimization of the simulation cell",&
         usage="TYPE (GEO_OPT|MD|DIRECT_CELL_OPT)",&
         enum_c_vals=s2a("GEO_OPT","MD","DIRECT_CELL_OPT"),&
         enum_desc=s2a(&
         "Performs a geometry optimization (the GEO_OPT section must be defined) between cell optimization steps. "//&
         " The stress tensor is computed at the optimized geometry.",&
         "Performs a molecular dynamics run (the MD section needs must defined) for computing the stress tensor "//&
         " used for the cell optimization.",&
         "Performs a geometry and cell optimization at the same time."//&
         " The stress tensor is computed at every step"),&
         enum_i_vals=(/default_cell_geo_opt_id,default_cell_md_id,default_cell_direct_id/),&
         default_i_val=default_cell_direct_id,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="EXTERNAL_PRESSURE",&
         description="Specifies the external pressure (1 value or the full 9 components of the pressure tensor) "//&
         "applied during the cell optimization.",&
         usage="EXTERNAL_PRESSURE {REAL} .. {REAL}",unit_str="bar",&
         default_r_vals=(/cp_unit_to_cp2k(100.0_dp,"bar",error=error),0.0_dp,0.0_dp,&
                          0.0_dp,cp_unit_to_cp2k(100.0_dp,"bar",error=error),0.0_dp,&
                          0.0_dp,0.0_dp,cp_unit_to_cp2k(100.0_dp,"bar",error=error)/),n_var=-1,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="KEEP_ANGLES",&
         description="Keep angles between the cell vectors constant, but allow the lenghts of the"//&
                     " cell vectors to change independently."//&
                     " Albeit general, this is most useful for triclinic cells, to enforce higher symmetry, see KEEP_SYMMETRY.",&
         usage="KEEP_ANGLES TRUE",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="KEEP_SYMMETRY",&
         description="Keep the requested initial cell symmetry (e.g. during a cell optimisation). "//&
                     "The initial symmetry must be specified in the &CELL section.",&
         usage="KEEP_SYMMETRY yes",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="PRESSURE_TOLERANCE",&
         description="Specifies the Pressure tolerance (compared to the external pressure) to achieve "//&
                     "during the cell optimization.",&
         usage="PRESSURE_TOLERANCE {REAL}",unit_str="bar",&
         default_r_val=cp_unit_to_cp2k(100.0_dp,"bar",error=error),error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    ! Create the PRINT subsection
    NULLIFY(subsection)
    CALL section_create(subsection,name="PRINT",&
         description="Controls the printing properties during a geometry optimization run",&
         n_keywords=0, n_subsections=1, repeats=.TRUE., required=.FALSE.,&
         error=error)
    NULLIFY(print_key)
    CALL cp_print_key_section_create(print_key,"program_run_info",&
         description="Controls the printing of basic information during the Geometry Optimization", &
         print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
         error=error)
    CALL section_add_subsection(subsection,print_key,error=error)
    CALL section_release(print_key,error=error)
    CALL cp_print_key_section_create(print_key,"cell",&
         description="Controls the printing of the cell eveytime a calculation using a new cell is started.", &
         print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
         unit_str="angstrom",error=error)
    CALL section_add_subsection(subsection,print_key,error=error)
    CALL section_release(print_key,error=error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

  END SUBROUTINE create_cell_opt_section

! *****************************************************************************
!> \brief creates the section for tuning transition states search
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] - University of Zurich - 01.2008
! *****************************************************************************
  SUBROUTINE create_ts_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ts_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection, &
                                                subsection2, subsection3

    failure=.FALSE.

    IF (.NOT.failure) THEN
       ! Create the Transition State subsection
       NULLIFY(section,keyword,subsection, subsection2)
       CALL section_create(section,name="TRANSITION_STATE",&
            description="Specifies parameters to perform a transition state search",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)

       CALL keyword_create(keyword, name="METHOD",&
            description="Specify which kind of method to use for locating transition states",&
            citations=(/Henkelman1999/),&
            usage="TYPE (DIMER)",&
            enum_c_vals=s2a("DIMER"),&
            enum_desc=s2a("Uses the dimer method to optimize transition states."),&
            enum_i_vals=(/default_dimer_method_id/),&
            default_i_val=default_dimer_method_id,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_create(subsection,name="DIMER",&
            description="Specifies parameters for Dimer Method",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)

       CALL keyword_create(keyword, name="DR",&
            description="This keyword sets the value for the DR parameter.",&
            usage="DR {real}",unit_str='angstrom',&
            default_r_val=cp_unit_to_cp2k(0.01_dp,"angstrom",error=error),error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="INTERPOLATE_GRADIENT",&
            description="This keyword controls the interpolation of the gradient whenever possible"//&
            " during the optimization of the Dimer. The use of this keywords saves 1 evaluation "//&
            " of energy/forces.", usage="INTERPOLATE_GRADIENT {logical}",default_l_val=.TRUE.,&
            lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ANGLE_TOLERANCE",&
            description="This keyword sets the value of the tolerance angle for the line search "//&
            " performed to optimize the orientation of the dimer.",&
            usage="ANGLE_TOL {real}",unit_str='rad',&
            default_r_val=cp_unit_to_cp2k(5.0_dp,"deg",error=error),error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_geoopt_section(subsection2, label="ROT_OPT",&
            description="This section sets the environment for the optimization of the rotation of the Dimer.",&
            just_optimizers=.TRUE.,error=error)
       NULLIFY(subsection3)
       CALL section_create(subsection3,name="PRINT",&
            description="Controls the printing properties during the dimer rotation optimization run",&
            n_keywords=0, n_subsections=1, repeats=.TRUE., required=.FALSE.,&
            error=error)
       NULLIFY(print_key)

       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="Controls the printing of basic information during the Geometry Optimization", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection3,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"ROTATIONAL_INFO",&
            description="Controls the printing basic info during the cleaning of the "//&
            "rotational degrees of freedom.", print_level=low_print_level,&
            add_last=add_last_numeric,filename="__STD_OUT__",error=error)
       CALL keyword_create(keyword, name="COORDINATES",&
            description="Prints atomic coordinates after rotation",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(subsection3,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL section_add_subsection(subsection2,subsection3,error=error)
       CALL section_release(subsection3,error=error)
       CALL section_add_subsection(subsection,subsection2,error=error)
       CALL section_release(subsection2,error=error)

       CALL section_create(subsection2,name="DIMER_VECTOR",&
            description="Specifies the initial dimer vector (used frequently to restart DIMER calculations)."//&
            " If not provided the starting orientation of the dimer is chosen randomly.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specify on each line the components of the dimer vector.",repeats=.TRUE.,&
            usage="{Real} {Real} {Real}", type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(subsection2,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(subsection,subsection2,error=error)
       CALL section_release(subsection2,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_ts_section

! *****************************************************************************
!> \brief creates the BFGS section
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] - University of Zurich - 01.2008
! *****************************************************************************
  SUBROUTINE create_bfgs_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_bfgs_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.
    IF (.NOT.failure) THEN
       ! create the BFGS subsection
       NULLIFY(section,keyword,print_key)
       CALL section_create(section,name="BFGS",&
            description="Provides parameters to tune the BFGS optimization",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)

       CALL keyword_create(keyword, name="TRUST_RADIUS",&
            description="Trust radius used in BFGS. Previously set to 0.1. "//&
            "Large values can lead to instabilities",&
            usage="TRUST_RADIUS {real}",unit_str='angstrom',&
            default_r_val=cp_unit_to_cp2k(0.25_dp,"angstrom",error=error),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="USE_MODEL_HESSIAN",&
            description="Uses a model Hessian as initial guess instead of a unit matrix."//&
            " Should lead in general to improved convergence might be switched off for exotic cases",&
            usage="USE_MODEL_HESSIAN",&
            default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="USE_RAT_FUN_OPT",&
            description="Includes a rational function optimization to determine the step."//&
            " Previously default but did not improve convergence in many cases",&
            usage="USE_RAT_FUN_OPT",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART_HESSIAN",&
            description="Controls the reading of the initial Hessian from file.",&
            usage="RESTART_HESSIAN",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART_FILE_NAME",&
            description="Specifies the name of the file used to read the initial Hessian.",&
            usage="RESTART_FILE_NAME {filename}",&
            default_lc_val="",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL cp_print_key_section_create(print_key,"RESTART",&
            description="Controls the printing of Hessian Restart file", &
            print_level=low_print_level,add_last=add_last_numeric,filename="BFGS",&
            common_iter_levels=2, error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_bfgs_section

! *****************************************************************************
!> \brief creates the CG section
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] - University of Zurich - 01.2008
! *****************************************************************************
  SUBROUTINE create_cg_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_cg_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection, subsubsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       ! create the CG subsection
       NULLIFY(section,subsection,subsubsection,keyword)
       CALL section_create(section,name="CG",&
            description="Provides parameters to tune the conjugate gradient optimization",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)

       CALL keyword_create(keyword, name="MAX_STEEP_STEPS",&
            description="Maximum number of steepest descent steps before starting the"//&
            " conjugate gradients optimization.",&
            usage="MAX_STEEP_STEPS {integer}",&
            default_i_val=0,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART_LIMIT",&
            description="Cosine of the angle between two consecutive searching directions."//&
            " If the angle during a CG optimization is less than the one corresponding to "//&
            " to the RESTART_LIMIT the CG is reset and one step of steepest descent is "//&
            " performed.",&
            usage="RESTART_LIMIT {real}",&
            default_r_val=0.9_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FLETCHER_REEVES",&
            description="Uses FLETCHER-REEVES instead of POLAK-RIBIERE when using Conjugate Gradients",&
            usage="FLETCHER-REEVES",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Line Search section
       CALL section_create(subsection,name="LINE_SEARCH",&
            description="Provides parameters to tune the line search during the conjugate gradient optimization",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)

       CALL keyword_create(keyword, name="TYPE",&
            description="1D line search algorithm to be used with the CG optimizer,"//&
                        " in increasing order of robustness and cost. ",&
            usage="TYPE GOLD",&
            default_i_val=ls_gold,&
            enum_c_vals=s2a( "NONE","2PNT","3PNT","GOLD","FIT"),&
            enum_desc=s2a("take fixed lenght steps",&
                          "extrapolate based on 2 points", &
                          "extrapolate based on on 3 points",&
                          "perform 1D golden section search of the minimum (very expensive)",&
                          "perform 1D fit of a parabola on several evaluation of energy "//&
                          "(very expensive and more robust vs numerical noise)"),&
            enum_i_vals=(/ls_none,ls_2pnt,ls_3pnt,ls_gold,ls_fit/),&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! 2PNT
       NULLIFY(subsubsection)
       CALL section_create(subsubsection,name="2PNT",&
            description="Provides parameters to tune the line search for the two point based line search.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)

       CALL keyword_create(keyword, name="MAX_ALLOWED_STEP",&
            description="Max allowed value for the line search step.",&
            usage="MAX_ALLOWED_STEP {real}",unit_str="internal_cp2k",&
            default_r_val=0.25_dp,error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LINMIN_GRAD_ONLY",&
            description="Use only the gradient, not the energy for line minimizations (e.g. in conjugate gradients).",&
            usage="LINMIN_GRAD_ONLY T",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(subsection,subsubsection,error=error)
       CALL section_release(subsubsection,error=error)

       ! GOLD
       NULLIFY(subsubsection)
       CALL section_create(subsubsection,name="GOLD",&
            description="Provides parameters to tune the line search for the gold search.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)

       CALL keyword_create(keyword, name="INITIAL_STEP",&
            description="Initial step size used, e.g. for bracketing or minimizers. "//&
            "Might need to be reduced for systems with close contacts",&
            usage="INITIAL_STEP {real}",unit_str="internal_cp2k",&
            default_r_val=0.2_dp,error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BRACK_LIMIT",&
            description="Limit in 1D bracketing during line search in Conjugate Gradients Optimization.",&
            usage="BRACK_LIMIT {real}",unit_str="internal_cp2k",&
            default_r_val=100.0_dp,error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BRENT_TOL",&
            description="Tolerance requested during Brent line search in Conjugate Gradients Optimization.",&
            usage="BRENT_TOL {real}",unit_str="internal_cp2k",&
            default_r_val=0.01_dp,error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BRENT_MAX_ITER",&
            description="Maximum number of iterations in brent algorithm "// &
            "(used for the line search in Conjugated Gradients Optimization)",&
            usage="BRENT_MAX_ITER {integer}",&
            default_i_val=100,error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(subsection,subsubsection,error=error)
       CALL section_release(subsubsection,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)
    END IF
  END SUBROUTINE create_cg_section

! *****************************************************************************
!> \brief creates the LBFGS section
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] - University of Zurich - 01.2008
! *****************************************************************************
  SUBROUTINE create_lbfgs_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_lbfgs_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    IF (.NOT.failure) THEN
       ! create the LBFGS section
       NULLIFY(section,keyword)
       CALL section_create(section,name="LBFGS",&
            description="Provides parameters to tune the limited memory BFGS (LBFGS) optimization",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            citations=(/Byrd1995/),&
            error=error)

       CALL keyword_create(keyword, name="MAX_H_RANK",&
            description="Maximum rank (and consequently size) of the "//&
            "approximate Hessian matrix used by the LBFGS optimizer. "//&
            "Larger values (e.g. 30) will accelerate the convergence behaviour "//&
            "at the cost of a larger memory consumption.",&
            usage="MAX_H_RANK {integer}",&
            default_i_val=5,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_F_PER_ITER",&
            description="Maximum number of force evaluations per iteration"// &
            "(used for the line search)",&
            usage="MAX_F_PER_ITER {integer}",&
            default_i_val=3,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WANTED_PROJ_GRADIENT",&
            description="Convergence criterium (overrides the general ones):"//&
            "Requested norm threshold of the gradient multiplied "// &
            "by the approximate Hessian.",&
            usage="WANTED_PROJ_GRADIENT {real}",unit_str="internal_cp2k",&
            default_r_val=1.0E-16_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WANTED_REL_F_ERROR",&
            description="Convergence criterium (overrides the general ones):"//&
            "Requested relative error on the objective function"//&
            "of the optimizer (the energy)",&
            usage="WANTED_REL_F_ERROR {real}",unit_str="internal_cp2k",&
            default_r_val=1.0E-16_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF

  END SUBROUTINE create_lbfgs_section

! *****************************************************************************
!> \brief creates the flexible_partitioning section
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Joost VandeVondele [04.2006]
! *****************************************************************************
  SUBROUTINE create_fp_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_fp_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="FLEXIBLE_PARTITIONING",&
            description="This section sets up flexible_partitioning",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword,print_key)

       CALL keyword_create(keyword, name="CENTRAL_ATOM",&
            description="Specifies the central atom.",&
            usage="CENTRAL_ATOM {integer}", required=.TRUE.,&
            n_var=1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="INNER_ATOMS",&
            description="Specifies the list of atoms that should remain close to the central atom.",&
            usage="INNER_ATOMS {integer} {integer} .. {integer}", required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="OUTER_ATOMS",&
            description="Specifies the list of atoms that should remain far from the central atom.",&
            usage="OUTER_ATOMS {integer} {integer} .. {integer}", required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="INNER_RADIUS",&
            description="radius of the inner wall",&
            usage="INNER_RADIUS {real} ", required=.TRUE., type_of_var=real_t, &
            n_var=1, unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="OUTER_RADIUS",&
            description="radius of the outer wall",&
            usage="OUTER_RADIUS {real} ", required=.TRUE., type_of_var=real_t, &
            n_var=1, unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="STRENGTH",&
            description="Sets the force constant of the repulsive harmonic potential",&
            usage="STRENGTH 1.0", default_r_val=1.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BIAS",&
            description="If a bias potential counter-acting the weight term should be applied (recommended).",&
            usage="BIAS F", default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TEMPERATURE",&
            description="Sets the temperature parameter that is used in the baising potential."//&
            "It is recommended to use the actual simulation temperature",&
            usage="TEMPERATURE 300", default_r_val=300.0_dp, unit_str='K', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SMOOTH_WIDTH",&
            description="Sets the width of the smooth counting function.",&
            usage="SMOOTH_WIDTH 0.2", default_r_val=0.02_dp, unit_str='angstrom', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL cp_print_key_section_create(print_key,"WEIGHTS",&
            description="Controls the printing of FP info during flexible partitioning simulations.", &
            print_level=low_print_level,common_iter_levels=1,&
            filename="FLEXIBLE_PARTIONING", error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"CONTROL",&
            description="Controls the printing of FP info at startup", &
            print_level=low_print_level,common_iter_levels=1, &
            filename="__STD_OUT__", error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF

  END SUBROUTINE create_fp_section


! *****************************************************************************
!> \brief ...
!> \param section will contain the driver section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author mceriotti
! *****************************************************************************
  SUBROUTINE create_driver_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_driver_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="DRIVER",&
            description="This section defines the parameters needed to run in i-PI driver mode.",&
            n_keywords=3, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword, name="unix",&
            description="Use a UNIX socket rather than an INET socket.",&
            usage="unix LOGICAL",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="port",&
            description="Port number for the i-PI server.",&
            usage="port <INTEGER>",&
            default_i_val=12345, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword, name="host",&
            description="Host name for the i-PI server.",&
            usage="host <HOSTNAME>",&
            default_c_val="localhost", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_driver_section

! *****************************************************************************
!> \brief creates the section for a path integral run
!> \param section will contain the pint section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE create_pint_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_pint_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection, &
                                                subsubsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="PINT",&
            description="The section that controls a path integral run",&
            n_keywords=11, n_subsections=8, repeats=.FALSE., required=.FALSE.,&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="p",&
            description="Specify number beads to use",repeats=.FALSE.,&
            default_i_val=3, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="proc_per_replica",&
            description="Specify number of processors to use for each replica",&
            repeats=.FALSE., default_i_val=0, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="num_steps",&
            description="Number of steps (if MAX_STEP is not explicitly given"//&
            " the program will perform this number of steps)",repeats=.FALSE.,&
            default_i_val=3, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="MAX_STEP",&
            description="Maximum step number (the program will stop if"//&
            " ITERATION >= MAX_STEP even if NUM_STEPS has not been reached)",&
            repeats=.FALSE., default_i_val=10, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="iteration",&
            description="Specify the iteration number from which it should be "//&
            "counted", default_i_val=0, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="Temp",&
            description="The temperature you want to simulate",&
            default_r_val=cp_unit_to_cp2k(300._dp,"K",error=error),&
            unit_str="K", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="T_tol",variants=(/"temp_to"/),&
            description="threshold for the oscillations of the temperature "//&
            "excedeed which the temperature is rescaled. 0 means no rescaling.",&
            default_r_val=0._dp,unit_str="K", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="dt",&
            description="timestep (might be subdivised in nrespa subtimesteps",&
            repeats=.FALSE.,&
            default_r_val=cp_unit_to_cp2k(1.0_dp,"fs",error=error),&
            usage="dt 1.0",unit_str="fs",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="nrespa",&
            description="number of respa steps for the bead for each md step",&
            repeats=.FALSE., default_i_val=5,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="transformation",&
            description="Specifies the coordinate transformation to use",&
            usage="TRANSFORMATION (NORMAL|STAGE)",&
            default_i_val=transformation_normal,&
            enum_c_vals=s2a("NORMAL","STAGE"),&
            enum_i_vals=(/transformation_normal,transformation_stage/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="FIX_CENTROID_POS",&
            description="Propagate all DOF but the centroid - "//&
            "useful for equilibration of the non-centroid modes "//&
            "(activated only if TRANSFORMATION==NORMAL)",&
            repeats=.FALSE., default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       NULLIFY(subsection,subsubsection)
       CALL section_create(subsection,name="NORMALMODE",&
            description="Controls the normal mode transformation",&
            n_keywords=3, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)
       CALL keyword_create(keyword, name="Q_CENTROID",&
            description="Value of the thermostat mass of centroid degree of freedom",&
            repeats=.FALSE., default_r_val=-1.0_dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="Q_BEAD",&
            description="Value of the thermostat mass of non-centroid degrees of freedom",&
            repeats=.FALSE., default_r_val=-1.0_dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="MODEFACTOR",&
            description="mass scale factor for non-centroid degrees of freedom",&
            repeats=.FALSE., default_r_val=1.0_dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="staging",&
            description="The section that controls the staging transformation",&
            n_keywords=2, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)
       CALL keyword_create(keyword, name="j",&
            description="Value of the j parameter for the staging transformation",&
            repeats=.FALSE., default_i_val=2,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="Q_END",&
            description="Value of the nose-hoover mass for the endbead (Q_end)",&
            repeats=.FALSE., default_i_val=2,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create( subsection, name="BEADS",&
            description="Sets positions and velocities of the beads",&
            n_keywords=0, n_subsections=2,&
            repeats=.FALSE., required=.FALSE., error=error )
       CALL create_coord_section(subsubsection,"BEADS",error=error)
       CALL section_add_subsection(subsection,subsubsection,error=error)
       CALL section_release(subsubsection,error=error)
       CALL create_velocity_section(subsubsection,"BEADS",error=error)
       CALL section_add_subsection(subsection,subsubsection,error=error)
       CALL section_release(subsubsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create( subsection, name="NOSE",&
            description="Controls the Nose-Hoover thermostats",&
            n_keywords=1, n_subsections=2,&
            repeats=.FALSE., required=.FALSE., error=error )
       CALL keyword_create(keyword, name="nnos",&
            description="length of nose-hoover chain. 0 means no thermostat",&
            repeats=.FALSE., default_i_val=2,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL create_coord_section(subsubsection,"NOSE",error=error)
       CALL section_add_subsection(subsection,subsubsection,error=error)
       CALL section_release(subsubsection,error=error)
       CALL create_velocity_section(subsubsection,"NOSE",error=error)
       CALL section_add_subsection(subsection,subsubsection,error=error)
       CALL section_release(subsubsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_gle_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL section_create( subsection, name="INIT",&
            description="Controls the initialization if the beads are not present",&
            repeats=.FALSE., required=.FALSE., error=error)

       CALL keyword_create(keyword, name="LEVY_POS_SAMPLE",&
            description="Sample bead positions assuming free particle "//&
            "behavior (performs a Levy random walk of length P around "//&
            "the classical position of each atom at the physical "//&
            "temperature defined in PINT%TEMP)",&
            repeats=.FALSE., default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="LEVY_CORRELATED",&
            description="Use the same Levy path for all atoms, though "//&
            "with mass-dependent variances (might help at very low T)",&
            repeats=.FALSE., default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="LEVY_TEMP_FACTOR",&
            description="Multiplicative correction factor for the "//&
            "temperature at which the Levy walk is performed "//&
            "(correction is due to the interactions that modify "//&
            "the spread of a free particle)",&
            repeats=.FALSE., default_r_val=1.0_dp,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,name="LEVY_SEED",&
            description="Initial seed for the (pseudo)random number "//&
            "generator that controls Levy walk for bead positions.",&
            usage="LEVY_SEED <INTEGER>",default_i_val=1234,&
            supported_feature=.TRUE.,repeats=.FALSE., error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="RANDOMIZE_POS",&
            description="add gaussian noise to the positions of the beads",&
            repeats=.FALSE., default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CENTROID_SPEED",&
            description="adds random velocity component to the centroid modes "//&
            "(useful to correct for the averaging out of the speed of various beads)",&
            repeats=.FALSE., default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="VELOCITY_QUENCH",&
            description="set the initial velocities to zero",&
            repeats=.FALSE., default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="VELOCITY_SCALE",&
            description="scale initial velocities to the temperature given in MOTION%PINT%TEMP",&
            repeats=.FALSE., default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_helium_section(subsection, error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="PRINT",&
            description="Controls the path integral-specific output",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(print_key)

       CALL cp_print_key_section_create(print_key,"ENERGY",&
            description="Controls the output of the path integral energies", &
            print_level=low_print_level,common_iter_levels=1,error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"CENTROID_POS",&
            description="Controls the output of the centroid's position", &
            unit_str="angstrom",&
            print_level=low_print_level,common_iter_levels=1,error=error)
       CALL add_format_keyword(keyword, print_key, pos=.TRUE.,&
            description="Output file format for the positions of centroid",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"CENTROID_VEL",&
            description="Controls the output of the centroid's velocity", &
            unit_str="bohr*au_t^-1",&
            print_level=low_print_level,common_iter_levels=1,error=error)
       CALL add_format_keyword(keyword, print_key, pos=.FALSE.,&
            description="Output file format for the velocity of centroid",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"CENTROID_GYR",&
            description="Controls the output of the centroid's radii of gyration", &
            unit_str="angstrom",&
            print_level=low_print_level,common_iter_levels=1,error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"COM",&
            description="Controls the output of the center of mass",&
            print_level=high_print_level,common_iter_levels=1,error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_pint_section

  ! ***************************************************************************
  !> \brief  Create the input section for superfluid helium solvent.
  !> \author Lukasz Walewski
  ! ***************************************************************************
! *****************************************************************************
!> \brief ...
!> \param section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_helium_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_helium_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection, &
                                                subsubsection

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (failure) THEN
      CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF

    CALL section_create(section,name="HELIUM",&
         description="The section that controls optional helium solvent"// &
         " environment (highly experimental, not for general use yet)",&
         n_keywords=11, n_subsections=4, repeats=.FALSE., required=.FALSE.,&
         error=error)

    NULLIFY(keyword)
    CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
         description="Whether or not to actually use this section",&
         usage="silent",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="HELIUM_ONLY",&
         description="Simulate helium solvent only, "//&
         "disregard solute entirely",&
         repeats=.FALSE., default_l_val=.FALSE.,&
         lone_keyword_l_val=.TRUE.,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="NUM_ENV", &
         description="Number of independent helium environments"// &
         " (only for restarts, do not set explicitly)",&
         repeats=.FALSE., default_i_val=1,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="POTENTIAL_FILE_NAME",&
         description="Name of the Helium interaction potential file",&
         repeats=.FALSE., default_lc_val="HELIUM.POT",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="NATOMS", &
         description="Number of helium atoms",&
         repeats=.FALSE., default_i_val=64,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="NBEADS", &
         description="Number of helium path integral beads",&
         repeats=.FALSE., default_i_val=25,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="INOROT", &
         description="Number of MC iterations at the same time slice(s)",&
         repeats=.FALSE., default_i_val=10000,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="IROT", &
         description="how often to reselect the time slice(s) to work on",&
         repeats=.FALSE., default_i_val=10000,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="BISECTION", &
         description="how many time slices to change at once (+1). "//&
         "Must be a power of 2 currently",&
         repeats=.FALSE., default_i_val=8,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MAX_PERM_CYCLE", &
         description="how large cyclic permutations to try",&
         repeats=.FALSE., default_i_val=6,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    NULLIFY(subsection)
    CALL section_create(subsection,name="M-SAMPLING",&
         description="Permutation cycle length sampling settings",&
         n_keywords=2, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
         error=error)
    CALL keyword_create(keyword, name="M-VALUE", &
         description="Value of m treated in a special way",&
         repeats=.FALSE.,&
         default_i_val=1,&
         error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="M-RATIO", &
         description="Probability ratio betw M-VALUE and other cycle lengths",&
         repeats=.FALSE.,&
         default_r_val=1.0_dp,&
         error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

    CALL keyword_create(keyword, name="PERIODIC", &
         description="Use periodic boundary conditions for helium",&
         repeats=.FALSE., default_l_val=.FALSE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="CELL_SIZE", &
         description="PBC unit cell size (NOTE 1: density, number of atoms"//&
         " and volume are interdependent - give only two of them; "//&
         "NOTE 2: for small cell sizes specify NATOMS instead)",&
         repeats=.FALSE.,type_of_var=real_t,unit_str="angstrom",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="CELL_SHAPE", &
         description="PBC unit cell shape for helium",&
         usage="CELL_SHAPE (CUBE|OCTAHEDRON)",&
         default_i_val=helium_cell_shape_cube,&
         enum_c_vals=s2a("CUBE","OCTAHEDRON"),&
         enum_i_vals=(/helium_cell_shape_cube,helium_cell_shape_octahedron/),&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="DENSITY", &
         description="trial density of helium for determining the helium "//&
         "box size",&
         repeats=.FALSE., &
         default_r_val=cp_unit_to_cp2k(0.02186_dp,"angstrom^-3",error=error),&
         unit_str="angstrom^-3",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="PRESAMPLE", &
         description="Presample He coordinates before first PIMD step",&
         repeats=.FALSE., default_l_val=.FALSE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="DROP_UNUSED_ENVS", &
         description="Drop He environments if N_restart > N_runtime "//&
         "(Warning: this will cause data loss in the restart file!)",&
         repeats=.FALSE., default_l_val=.FALSE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL section_create(subsection,name="RDF",&
         description="Radial distribution function generation settings",&
         n_keywords=2, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
         error=error)
    CALL keyword_create(keyword, name="MAXR", &
         description="Maximum RDF range, defaults to unit cell size",&
         repeats=.FALSE.,type_of_var=real_t,&
         unit_str="angstrom",error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="NBIN", &
         description="Number of bins",&
         repeats=.FALSE.,&
         default_i_val=700,&
         error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

    NULLIFY(subsection)
    CALL section_create(subsection,name="RHO",&
         description="Density distribution settings",&
         n_keywords=2, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
         error=error)
    CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
         description="Whether or not to actually calculate densities "//&
         "(requires significant amount of memory, depending on the value of NBIN)",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,&
         error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="NBIN", &
         description="Number of bins",&
         repeats=.FALSE.,&
         default_i_val=100,&
         error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="IWEIGHT", &
         description="Weight the restarted density should be given " // &
         "(number of MC steps used to average the restarted density, " // &
         "negative value - the same weight as the run-time density, " // &
         "usually should not be changed)",&
         repeats=.FALSE.,&
         default_i_val=-1,&
         error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    NULLIFY(subsubsection)
    CALL section_create(subsubsection,name="CUBE_DATA",&
         description="Density data used for restarts",&
         n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
         error=error)
    CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
         description="Cubefile data", &
         repeats=.TRUE., usage="{Real} ...",&
         type_of_var=real_t, n_var=-1, error=error)
    CALL section_add_keyword(subsubsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL section_add_subsection(subsection,subsubsection,error=error)
    CALL section_release(subsubsection,error=error)

    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)
    ! end of subsection RHO

    CALL create_coord_section(subsection,"HELIUM",error=error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

    CALL section_create(subsection,name="PERM",&
         description="Permutation state used for restart",&
         n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
         error=error)
    CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
         description="Specify particle index permutation for every "// &
         "helium atom",repeats=.TRUE.,usage="{Integer} ...",&
         type_of_var=integer_t, n_var=-1, error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

    CALL section_create(subsection,name="FORCE",&
         description="Forces exerted by the helium on the solute system"//&
         " (used for restarts)",&
         n_keywords=0, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)
    CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_", &
         description="Number of real values should be 3 * "//&
         "<num_solute_atoms> * <num_solute_beads>", repeats=.TRUE., &
         usage="{Real} ...", type_of_var=real_t, &
         n_var=-1, error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

    CALL section_create(subsection,name="RNG_STATE",&
         description="Random number generator state for all processors",&
         n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
         error=error)
    CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
         description="Three real arrays of DIMENSION(3,2) times two RNG "//&
         "streams - 36 real values per processor",&
         repeats=.TRUE.,usage="automatically filled, do not edit by hand",&
         type_of_var=real_t, n_var=-1, error=error)
    CALL section_add_keyword(subsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

    CALL section_create(subsection,name="PRINT",&
         description="The section that controls the output of the helium code",&
         n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
         error=error)

    ! *************************************************************************
    !> Printkeys for properites output
    ! *************************************************************************
    NULLIFY(print_key)

    ! Properties printed at SILENT print level
    !

    ! Properties printed at LOW print level
    !
    CALL cp_print_key_section_create(print_key,"ENERGY",&
         description="Controls the output of the helium energies"//&
         " (averaged over MC step)", &
         print_level=low_print_level,common_iter_levels=1,error=error)
    CALL section_add_subsection(subsection,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"SDENSITY",&
         description="Controls the output of the helium superfluid density",&
         print_level=low_print_level,common_iter_levels=1,error=error)
    CALL section_add_subsection(subsection,print_key,error=error)
    CALL section_release(print_key,error=error)

    ! Properties printed at MEDIUM print level
    !
    CALL cp_print_key_section_create(print_key,"COORDINATES",&
         description="Controls the output of helium coordinates",&
         print_level=medium_print_level,common_iter_levels=1,error=error)
    CALL keyword_create(keyword, name="FORMAT",&
         description="Output file format for the coordinates",&
         usage="FORMAT (PDB|XYZ)",&
         default_i_val=fmt_id_pdb,&
         enum_c_vals=s2a("PDB","XYZ"),&
         enum_i_vals= (/fmt_id_pdb,fmt_id_xyz/),&
         enum_desc=s2a( "Bead coordinates and connectivity is written in PDB format",&
                        "Only bead coordinates are written in XYZ format"),&
            error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL section_add_subsection(subsection,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"RDF",&
         description="Controls the output of the helium radial distribution function",&
         print_level=medium_print_level,common_iter_levels=1,error=error)
    CALL section_add_subsection(subsection,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"RHO",&
         description="Controls the output of the helium density "//&
                     "(Gaussian cube file format)",&
         each_iter_names=s2a("MD"),each_iter_values=(/100/),&
         print_level=medium_print_level, common_iter_levels=1,&
         add_last=add_last_numeric, error=error)
    CALL keyword_create(keyword, name="BACKUP_COPIES",&
         description="Specifies the maximum number of backup copies.",&
         usage="BACKUP_COPIES {int}",&
         default_i_val=1, error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL section_add_subsection(subsection,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"PLENGTH",&
         description="Controls the output of the helium permutation length",&
         print_level=medium_print_level,common_iter_levels=1,error=error)
    CALL section_add_subsection(subsection,print_key,error=error)
    CALL section_release(print_key,error=error)

    ! Properties printed at HIGH print level
    !
    CALL cp_print_key_section_create(print_key,"ACCEPTS",&
         description="Controls the output of the helium acceptance data",&
         print_level=high_print_level,common_iter_levels=1,error=error)
    CALL section_add_subsection(subsection,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"PERM",&
         description="Controls the output of the helium permutation state",&
         print_level=high_print_level,common_iter_levels=1,error=error)
    CALL section_add_subsection(subsection,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"FORCES",&
         description="Controls the output of the helium forces on the solute",&
         print_level=high_print_level,common_iter_levels=1,error=error)
    CALL section_add_subsection(subsection,print_key,error=error)
    CALL section_release(print_key,error=error)

    ! Properties printed at DEBUG print level
    !
    CALL cp_print_key_section_create(print_key,"FORCES_INST",&
         description="Controls the output of the instantaneous helium forces on the solute",&
         print_level=debug_print_level,common_iter_levels=1,error=error)
    CALL section_add_subsection(subsection,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"WNUMBER",&
         description="Controls the output of the helium winding number",&
         print_level=debug_print_level,common_iter_levels=1,error=error)
    CALL section_add_subsection(subsection,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

    RETURN
  END SUBROUTINE create_helium_section


END MODULE input_cp2k_motion
