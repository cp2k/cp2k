!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!


! *****************************************************************************
!> \par History
!>      - taken out of input_cp2k_motion
!> \author Ole Schuett
! *****************************************************************************

MODULE input_cp2k_neb
  USE bibliography,                    ONLY: Elber1987,&
                                             Jonsson1998,&
                                             Jonsson2000_1,&
                                             Jonsson2000_2,&
                                             Wales2004
  USE cp_output_handling,              ONLY: add_last_numeric,&
                                             cp_print_key_section_create,&
                                             high_print_level,&
                                             low_print_level,&
                                             medium_print_level
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE input_constants,                 ONLY: &
       band_diis_opt, band_md_opt, do_b_neb, do_ci_neb, do_d_neb, do_eb, &
       do_it_neb, do_rep_blocked, do_rep_interleaved, do_sm, pot_neb_fe, &
       pot_neb_full, pot_neb_me
  USE input_cp2k_thermostats,          ONLY: create_coord_section,&
                                             create_velocity_section
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_neb'

PUBLIC :: create_band_section


CONTAINS

! *****************************************************************************
!> \brief creates the section for a BAND run
!> \param section will contain the pint section
!> \author Teodoro Laino 09.2006 [tlaino]
! *****************************************************************************
  SUBROUTINE create_band_section(section)
    TYPE(section_type), POINTER              :: section

    CHARACTER(len=*), PARAMETER :: routineN = 'create_band_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection, &
                                                subsubsection

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,failure)
    CALL section_create(section,name="band",&
         description="The section that controls a BAND run",&
         n_keywords=1, n_subsections=0, repeats=.FALSE., &
         citations=(/Elber1987,Jonsson1998,Jonsson2000_1,Jonsson2000_2,Wales2004/))
    NULLIFY(keyword, print_key, subsection, subsubsection)

    CALL keyword_create(keyword, name="NPROC_REP",&
         description="Specify the number of processors to be used per replica "//&
         "environment (for parallel runs)",&
         default_i_val=1)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="PROC_DIST_TYPE",&
         description="Specify the topology of the mapping of processors into replicas.",&
         usage="PROC_DIST_TYPE (INTERLEAVED|BLOCKED)",&
         enum_c_vals=s2a("INTERLEAVED",&
                         "BLOCKED"),&
         enum_desc=s2a( "Interleaved distribution",&
                        "Blocked distribution"),&
         enum_i_vals=(/do_rep_interleaved,do_rep_blocked/),&
         default_i_val=do_rep_blocked)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="BAND_TYPE",&
         description="Specifies the type of BAND calculation",&
         usage="BAND_TYPE (B-NEB|IT-NEB|CI-NEB|D-NEB|SM|EB)",&
         default_i_val=do_it_neb,&
         enum_c_vals=s2a("B-NEB",&
                         "IT-NEB",&
                         "CI-NEB",&
                         "D-NEB",&
                         "SM",&
                         "EB"),&
         enum_desc=s2a( "Bisection nudged elastic band",&
                        "Improved tangent nudged elastic band",&
                        "Climbing image nudged elastic band",&
                        "Doubly nudged elastic band",&
                        "String Method",&
                        "Elastic band (Hamiltonian formulation)"),&
         enum_i_vals=(/do_b_neb,do_it_neb,do_ci_neb,do_d_neb,do_sm,do_eb/))
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="NUMBER_OF_REPLICA",&
         description="Specify the number of Replica to use in the BAND",&
         default_i_val=10)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="USE_COLVARS",&
         description="Uses a version of the band scheme projected in a subspace of colvars.",&
         default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="POT_TYPE",&
         description="Specifies the type of potential used in the BAND calculation",&
         usage="POT_TYPE (FULL|FE|ME)",&
         default_i_val=pot_neb_full,&
         enum_c_vals=s2a("FULL",&
                         "FE",&
                         "ME"),&
         enum_desc=s2a( "Full potential (no projections in a subspace of colvars)",&
                        "Free energy (requires a projections in a subspace of colvars)",&
                        "Minimum energy (requires a projections in a subspace of colvars)"),&
         enum_i_vals=(/pot_neb_full,pot_neb_fe,pot_neb_me/))
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="ROTATE_FRAMES",&
         description="Compute at each BAND step the RMSD and rotate the frames in order"//&
         " to minimize it.",&
         default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="ALIGN_FRAMES",&
         description="Enables the alignment of the frames at the beginning of a BAND calculation. "//&
         "This keyword does not affect the rotation of the replicas during a BAND calculation.",&
         default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="K_SPRING",&
         variants=(/"K"/),&
         description="Specify the value of the spring constant",&
         default_r_val=0.02_dp)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    ! Convergence_control
    CALL section_create(subsection,name="CONVERGENCE_CONTROL",&
         description="Setup parameters to control the convergence criteria for BAND",&
         repeats=.FALSE.)
    CALL keyword_create(keyword, name="MAX_DR",&
         description="Tolerance on the maximum value of the displacement on the BAND.",&
         usage="MAX_DR {real}",&
         default_r_val=0.0002_dp)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="MAX_FORCE",&
         description="Tolerance on the maximum value of Forces on the BAND.",&
         usage="MAX_FORCE {real}",&
         default_r_val=0.00045_dp)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="RMS_DR",&
         description="Tolerance on RMS displacements on the BAND.",&
         usage="RMS_DR {real}",&
         default_r_val=0.0001_dp)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="RMS_FORCE",&
         description="Tolerance on RMS Forces on the BAND.",&
         usage="RMS_FORCE {real}",&
         default_r_val=0.00030_dp)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)
    CALL section_add_subsection(section, subsection)
    CALL section_release(subsection)

    NULLIFY(subsection, subsubsection)
    ! CI-NEB section
    CALL section_create(subsection,name="CI_NEB",&
         description="Controls parameters for CI-NEB type calculation only.",&
         repeats=.FALSE.)
    CALL keyword_create(keyword, name="NSTEPS_IT",&
         description="Specify the number of steps of IT-NEB to perform before "//&
         "switching on the CI algorithm",&
         default_i_val=5)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)
    CALL section_add_subsection(section, subsection)
    CALL section_release(subsection)

    ! String Method section
    CALL section_create(subsection,name="STRING_METHOD",&
         description="Controls parameters for String Method type calculation only.",&
         repeats=.FALSE.)

    CALL keyword_create(keyword, name="SPLINE_ORDER",&
         description="Specify the oder of the spline used in the String Method.",&
         default_i_val=1)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)
    CALL keyword_create(keyword, name="SMOOTHING",&
         description="Smoothing parameter for the reparametrization of the frames.",&
         default_r_val=0.2_dp)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL section_add_subsection(section, subsection)
    CALL section_release(subsection)

    ! Optimization section
    CALL section_create(subsection,name="optimize_band",&
         description="Specify the optimization method for the band",&
         repeats=.TRUE.)
    CALL create_opt_band_section(subsection)
    CALL section_add_subsection(section, subsection)
    CALL section_release(subsection)

    ! replica section: to specify coordinates and velocities (possibly) of the
    ! different replica used in the BAND
    CALL section_create(subsection,name="replica",&
         description="Specify coordinates and velocities (possibly) of the replica",&
         repeats=.TRUE.)
    ! Colvar
    CALL keyword_create(keyword, name="COLLECTIVE",&
         description="Specifies the value of the collective variables used in the projected"//&
         " BAND method. The order of the values is the order of the COLLECTIVE section in the"//&
         " constraints/restraints section",&
         usage="COLLECTIVE {real} .. {real}",&
         type_of_var=real_t, n_var=-1)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)
    ! Coordinates read through an external file
    CALL keyword_create(keyword, name="COORD_FILE_NAME",&
         description="Name of the xyz file with coordinates (alternative to &COORD section)",&
         usage="COORD_FILE_NAME <CHAR>",&
         default_lc_val="")
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)
    ! Coordinates and velocities
    CALL create_coord_section(subsubsection,"BAND")
    CALL section_add_subsection(subsection,subsubsection)
    CALL section_release(subsubsection)
    CALL create_velocity_section(subsubsection,"BAND")
    CALL section_add_subsection(subsection,subsubsection)
    CALL section_release(subsubsection)

    CALL section_add_subsection(section, subsection)
    CALL section_release(subsection)

    ! Print key section
    CALL cp_print_key_section_create(print_key,"program_run_info",&
         description="Controls the printing basic info about the BAND run", &
         print_level=medium_print_level,add_last=add_last_numeric,filename="__STD_OUT__")

    CALL keyword_create(keyword, name="INITIAL_CONFIGURATION_INFO",&
         description="Print information for the setup of the initial configuration.",&
         usage="INITIAL_CONFIGURATION_INFO <LOGICAL>",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(print_key,keyword)
    CALL keyword_release(keyword)

    CALL section_add_subsection(section,print_key)
    CALL section_release(print_key)

    CALL cp_print_key_section_create(print_key,"convergence_info",&
         description="Controls the printing of the convergence criteria during a BAND run", &
         print_level=medium_print_level,add_last=add_last_numeric,filename="__STD_OUT__")
    CALL section_add_subsection(section,print_key)
    CALL section_release(print_key)

    CALL cp_print_key_section_create(print_key,"replica_info",&
         description="Controls the printing of each replica info during a BAND run", &
         print_level=medium_print_level,add_last=add_last_numeric,filename="__STD_OUT__")
    CALL section_add_subsection(section,print_key)
    CALL section_release(print_key)

    CALL cp_print_key_section_create(print_key,"ENERGY",&
         description="Controls the printing of the ENER file in a BAND run", &
         print_level=low_print_level, common_iter_levels=1,&
         filename="")
    CALL section_add_subsection(section,print_key)
    CALL section_release(print_key)

    CALL cp_print_key_section_create(print_key,"BANNER",&
         description="Controls the printing of the BAND banner", &
         print_level=low_print_level, common_iter_levels=1,&
         filename="__STD_OUT__")
    CALL section_add_subsection(section,print_key)
    CALL section_release(print_key)
  END SUBROUTINE create_band_section

! *****************************************************************************
!> \brief creates the optimization section for a BAND run
!> \param section will contain the pint section
!> \author Teodoro Laino 02.2007 [tlaino]
! *****************************************************************************
  SUBROUTINE create_opt_band_section(section)
    TYPE(section_type), POINTER              :: section

    CHARACTER(len=*), PARAMETER :: routineN = 'create_opt_band_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection, &
                                                subsubsection

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,failure)
    NULLIFY(keyword, print_key, subsection, subsubsection)

    CALL keyword_create(keyword, name="OPT_TYPE",&
         description="Specifies the type optimizer used for the band",&
         usage="OPT_TYPE (MD|DIIS)",&
         default_i_val=band_diis_opt,&
         enum_c_vals=s2a("MD",&
                         "DIIS"),&
         enum_desc=s2a( "Molecular dynamics-based optimizer",&
                        "Coupled steepest descent / direct inversion in the iterative subspace"),&
         enum_i_vals=(/band_md_opt,band_diis_opt/))
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="OPTIMIZE_END_POINTS",&
         description="Performs also an optimization of the end points of the band.",&
         default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(section,keyword)
    CALL keyword_release(keyword)

    ! MD optimization section
    CALL section_create(subsection,name="MD",&
         description="Activate the MD based optimization procedure for BAND",&
         repeats=.FALSE.)

    CALL keyword_create(keyword, name="MAX_STEPS",&
         description="Specify the maximum number of MD steps",&
         default_i_val=100)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(&
         keyword=keyword,&
         name="timestep",&
         description="The length of an integration step",&
         usage="timestep 1.0",&
         default_r_val=cp_unit_to_cp2k(value=0.5_dp,&
         unit_str="fs"),&
         unit_str="fs")
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="TEMPERATURE",&
         description="Specify the initial temperature",&
         default_r_val=cp_unit_to_cp2k(value=0.0_dp,&
         unit_str="K"),&
         unit_str="K")
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    ! Temp_control
    CALL section_create(subsubsection,name="TEMP_CONTROL",&
         description="Setup parameters to control the temperature during a BAND MD run.",&
         repeats=.FALSE.)
    CALL keyword_create(keyword, name="TEMPERATURE",&
         description="Specify the target temperature",&
         type_of_var=real_t,unit_str="K")
    CALL section_add_keyword(subsubsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="TEMP_TOL",&
         description="Specify the tolerance on the temperature for rescaling",&
         default_r_val=cp_unit_to_cp2k(value=0.0_dp,&
         unit_str="K"),&
         unit_str="K")
    CALL section_add_keyword(subsubsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="TEMP_TOL_STEPS",&
         description="Specify the number of steps to apply a temperature control",&
         default_i_val=0)
    CALL section_add_keyword(subsubsection,keyword)
    CALL keyword_release(keyword)
    CALL section_add_subsection(subsection, subsubsection)
    CALL section_release(subsubsection)

    ! Vel_control
    CALL section_create(subsubsection,name="VEL_CONTROL",&
         description="Setup parameters to control the velocity during a BAND MD run.",&
         repeats=.FALSE.)
    CALL keyword_create(keyword, name="ANNEALING",&
         description="Specify the annealing coefficient",&
         default_r_val=1.0_dp)
    CALL section_add_keyword(subsubsection,keyword)
    CALL keyword_release(keyword)
    CALL keyword_create(keyword, name="PROJ_VELOCITY_VERLET",&
         description="Uses a Projected Velocity Verlet instead of a normal Velocity Verlet."//&
         " Every time the cosine between velocities and forces is < 0 velocities are"//&
         " zeroed.",&
         usage="PROJ_VELOCITY_VERLET <LOGICAL>",&
         default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(subsubsection,keyword)
    CALL keyword_release(keyword)
    CALL keyword_create(keyword, name="SD_LIKE",&
         description="Zeros velocity at each MD step emulating a steepest descent like"//&
         "(SD_LIKE) approach",&
         usage="SD_LIKE <LOGICAL>",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(subsubsection,keyword)
    CALL keyword_release(keyword)
    CALL section_add_subsection(subsection, subsubsection)
    CALL section_release(subsubsection)
    ! End of MD
    CALL section_add_subsection(section, subsection)
    CALL section_release(subsection)

    ! DIIS optimization section
    CALL section_create(subsection,name="DIIS",&
         description="Activate the DIIS based optimization procedure for BAND",&
         repeats=.FALSE.)

    CALL keyword_create(keyword, name="MAX_SD_STEPS",&
         description="Specify the maximum number of SD steps to perform"//&
         " before switching on DIIS (the minimum number will always be equal to N_DIIS).",&
         default_i_val=1)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="MAX_STEPS",&
         description="Specify the maximum number of optimization steps",&
         default_i_val=100)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="N_DIIS",&
         variants=(/"NDIIS"/),&
         description="Number of history vectors to be used with DIIS",&
         usage="N_DIIS 4",&
         default_i_val=7)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="STEPSIZE",&
         description="Initial stepsize used for the line search, sometimes this parameter"//&
         "can be reduced to stablize DIIS",&
         usage="STEPSIZE <REAL>",&
         default_r_val=1.0_dp)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="MAX_STEPSIZE",&
         description="Maximum stepsize used for the line search, sometimes this parameter"//&
         "can be reduced to stablize the LS for particularly difficult initial geometries",&
         usage="MAX_STEPSIZE <REAL>",&
         default_r_val=2.0_dp)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="NP_LS",&
         description="Number of points used in the line search SD.",&
         usage="NP_LS <INTEGER>",&
         default_i_val=2)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="NO_LS",&
         description="Does not perform LS during SD. Useful in combination with a proper STEPSIZE"//&
         " for particularly out of equilibrium starting geometries.",&
         default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name="CHECK_DIIS",&
         description="Performes a series of checks on the DIIS solution in order to accept the DIIS step."//&
         " If set to .FALSE. the only check performed is that the angle between the DIIS solution and the"//&
         " reference vector is less than Pi/2. Can be useful if many DIIS steps are rejected.",&
         default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
    CALL section_add_keyword(subsection,keyword)
    CALL keyword_release(keyword)

    CALL cp_print_key_section_create(print_key,"diis_info",&
         description="Controls the printing of DIIS info during a BAND run", &
         print_level=high_print_level,add_last=add_last_numeric,filename="__STD_OUT__")
    CALL section_add_subsection(subsection,print_key)
    CALL section_release(print_key)

    CALL section_add_subsection(section, subsection)
    CALL section_release(subsection)
  END SUBROUTINE create_opt_band_section


END MODULE input_cp2k_neb
