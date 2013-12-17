!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      10.2005 split input_cp2k into smaller modules [fawzi]
!> \author teo & fawzi
! *****************************************************************************
MODULE input_cp2k_motion
  USE bibliography,                    ONLY: &
       Byrd1995, Elber1987, Evans1983, Guidon2008, Henkelman1999, &
       Jonsson1998, Jonsson2000_1, Jonsson2000_2, Kolafa2004, Kuhne2007, &
       Minary2003, Ricci2003, Tuckerman1992, VandeVondele2002, Wales2004
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE f77_blas
  USE input_constants
  USE input_cp2k_barostats,            ONLY: create_barostat_section
  USE input_cp2k_constraints,          ONLY: create_constraint_section
  USE input_cp2k_free_energy,          ONLY: create_fe_section
  USE input_cp2k_subsys,               ONLY: create_structure_data_section
  USE input_cp2k_thermostats,          ONLY: create_coord_section,&
                                             create_gle_section,&
                                             create_region_section,&
                                             create_thermo_fast_section,&
                                             create_thermo_slow_section,&
                                             create_thermostat_section,&
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
                                             lchar_t,&
                                             real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_motion'

PUBLIC :: create_motion_section,&
          create_band_section,&
          add_format_keyword,&
          create_md_section

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
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="motion",&
            description="This section defines a set of tool connected with the motion of the nuclei.",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection, print_key)

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

       CALL section_create(subsection,name="print",&
            description="Controls the printing properties during an MD run",&
            n_keywords=0, n_subsections=1, repeats=.TRUE., required=.FALSE.,&
            error=error)
       NULLIFY(keyword)
       CALL cp_print_key_section_create(print_key,"TRAJECTORY",&
            description="Controls the output of the trajectory",&
            print_level=low_print_level, common_iter_levels=1,&
            filename="",unit_str="angstrom",error=error)
       CALL add_format_keyword(keyword, print_key, pos=.TRUE.,&
            description="Specifies the format of the output file for the trajectory.", error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SHELL_TRAJECTORY",&
            description="Controls the output of the trajectory of shells when the shell-model is used ",&
            print_level=high_print_level, common_iter_levels=1,&
            filename="",unit_str="angstrom",error=error)
       CALL add_format_keyword(keyword, print_key, pos=.TRUE.,&
            description="Specifies the format of the output file for the trajectory of shells.", error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"CORE_TRAJECTORY",&
            description="Controls the output of the trajectory of cores when the shell-model is used ",&
            print_level=high_print_level, common_iter_levels=1,&
            filename="",unit_str="angstrom",error=error)
       CALL add_format_keyword(keyword, print_key, pos=.TRUE.,&
            description="Specifies the format of the output file for the trajectory of cores.", error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"CELL",&
            description="Controls the output of the simulation cell. "//&
            "For later analysis of the trajectory it is recommendable that the "//&
            "frequency of printing is the same as the one used for the trajectory file.",&
            print_level=high_print_level, common_iter_levels=1,&
            filename="",error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"VELOCITIES",&
            description="Controls the output of the velocities",&
            print_level=high_print_level, common_iter_levels=1,&
            filename="",unit_str="bohr*au_t^-1",error=error)
       CALL add_format_keyword(keyword, print_key, pos=.FALSE.,&
            description="Specifies the format of the output file for the velocities.", error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SHELL_VELOCITIES",&
            description="Controls the output of the velocities of shells when the shell model is used",&
            print_level=high_print_level, common_iter_levels=1,&
            filename="",unit_str="bohr*au_t^-1",error=error)
       CALL add_format_keyword(keyword, print_key, pos=.FALSE.,&
            description="Specifies the format of the output file for the velocities of shells.", error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"CORE_VELOCITIES",&
            description="controls the output of the velocities of cores when the shell model is used",&
            print_level=high_print_level, common_iter_levels=1,&
            filename="",unit_str="bohr*au_t^-1",error=error)
       CALL add_format_keyword(keyword, print_key, pos=.FALSE.,&
            description="Specifies the format of the output file for the velocities of cores.", error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL create_structure_data_section(print_key, error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"FORCE_MIXING_LABELS",&
            description="Controls the output of the force mixing (FORCE_EVAL&QMMM&FORCE_MIXING) labels",&
            print_level=high_print_level, common_iter_levels=1,&
            filename="",error=error)
       CALL add_format_keyword(keyword, print_key, pos=.FALSE.,&
            description="Specifies the format of the output file for the force mixing labels.", error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"FORCES",&
            description="Controls the output of the forces",&
            print_level=high_print_level, common_iter_levels=1,&
            filename="",unit_str="hartree*bohr^-1",error=error)
       CALL add_format_keyword(keyword, print_key, pos=.FALSE.,&
            description="Specifies the format of the output file for the forces.", error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SHELL_FORCES",&
            description="controls the output of the forces on shells when shell-model is used",&
            print_level=high_print_level, common_iter_levels=1,&
            filename="",unit_str="hartree*bohr^-1",error=error)
       CALL add_format_keyword(keyword, print_key, pos=.FALSE.,&
            description="Specifies the format of the output file for the forces on shells.", error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"CORE_FORCES",&
            description="controls the output of the forces on cores when shell-model is used",&
            print_level=high_print_level, common_iter_levels=1,&
            filename="",unit_str="hartree*bohr^-1",error=error)
       CALL add_format_keyword(keyword, print_key, pos=.FALSE.,&
            description="Specifies the format of the output file for the forces on cores.", error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"MIXED_ENERGIES",&
            description="Controls the output of the energies of the two"//&
            "regular FORCE_EVALS in the MIXED method"//&
            "printed is step,time,Etot,E_F1,E_F2,CONS_QNT",&
            print_level=low_print_level, common_iter_levels=1,&
            filename="",error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"STRESS",&
            description="Controls the output of the stress tensor",&
            print_level=high_print_level, common_iter_levels=1,&
            filename="",error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"RESTART",&
            description="Controls the dumping of the restart file during runs. "//&
            "By default keeps a short history of three restarts. See also RESTART_HISTORY", &
            each_iter_names=s2a("MD"),each_iter_values=(/20/), &
            print_level=silent_print_level, common_iter_levels=1,  &
            add_last=add_last_numeric, filename="",error=error)

       CALL keyword_create(keyword, name="BACKUP_COPIES",&
            description="Specifies the maximum number of backup copies.",&
            usage="BACKUP_COPIES {int}",&
            default_i_val=3, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SPLIT_RESTART_FILE",&
                           description="If specified selected input sections, which are growing with the "//&
                                       "number of atoms in the system, are written to another restart file "//&
                                       "in binary format instead of the default restart file in human "//&
                                       "readable ASCII format. This split of the restart file may "//&
                                       "provide significant memory savings and an accelerated I/O for "//&
                                       "systems with a very large number of atoms",&
                           usage="SPLIT_RESTART_FILE yes",&
                           default_l_val=.FALSE.,&
                           lone_keyword_l_val=.TRUE.,&
                           error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"RESTART_HISTORY",&
            description="Dumps unique restart files during the run keeping all of them."//&
            "Most useful if recovery is needed at a later point.",&
            print_level=low_print_level, common_iter_levels=0, &
            each_iter_names=s2a("MD","GEO_OPT","ROT_OPT"),each_iter_values=(/500,500,500/), &
            filename="",error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)


       CALL cp_print_key_section_create(print_key,"TRANSLATION_VECTOR",&
            description="Dumps the translation vector applied along an MD (if any). Useful"//&
            " for postprocessing of QMMM trajectories in which the QM fragment is continuously"//&
            " centered in the QM box",&
            print_level=high_print_level, common_iter_levels=1,&
            filename="",error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

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
                     " For example, a cubic cell might become orthorhombic, but not triclinic.",&
         usage="KEEP_ANGLES TRUE",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="KEEP_SYMMETRY",&
         description="Keep the requested initial cell symmetry (e.g. during a cell optimisation). "//&
                     "See cell section for how to impose a cell symmetry.",&
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
            usage="<REAL> <REAL> <REAL>", type_of_var=real_t, n_var=-1, error=error)
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
!> \param section will contain the md section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE create_md_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_md_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="MD",&
            description="This section defines the whole set of parameters needed perform an MD run.",&
            n_keywords=13, n_subsections=6, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword, subsection)
       CALL keyword_create(keyword, name="ensemble",&
            description="The ensemble/integrator that you want to use for MD propagation",&
            usage="ensemble nve",&
            default_i_val=nve_ensemble,&
            enum_c_vals=s2a("NVE", "NVT", "NPT_I", "NPT_F", "MSST", "MSST_DAMPED",&
            "HYDROSTATICSHOCK", "ISOKIN", "REFTRAJ", "LANGEVIN", "NPE_F", &
            "NPE_I", "NVT_ADIABATIC" ),&
            enum_desc=s2a("constant energy (microcanonical)",&
            "constant temperature and volume (canonical)",&
            "constant temperature and pressure using an isotropic cell",&
            "constant temperature and pressure using a flexible cell",&
            "simulate steady shock (uniaxial)",&
            "simulate steady shock (uniaxial) with extra viscosity",&
            "simulate steady shock with hydrostatic pressure",&
            "constant kinetic energy",&
            "reading frames from a file called reftraj.xyz (e.g. for property calculation)",&
            "langevin dynamics (constant temperature)",&
            "constant pressure ensemble (no thermostat)",&
            "constant pressure ensemble using an isotropic cell (no thermostat)",&
            "adiabatic dynamics in constant temperature and volume ensemble (CAFES)"),&
            citations=(/Evans1983,VandeVondele2002,Minary2003,Kolafa2004/),&
            enum_i_vals=(/nve_ensemble,nvt_ensemble,npt_i_ensemble,npt_f_ensemble,&
            nph_uniaxial_ensemble,nph_uniaxial_damped_ensemble,nph_ensemble,isokin_ensemble,&
            reftraj_ensemble,langevin_ensemble, npe_f_ensemble, npe_i_ensemble, &
            nvt_adiabatic_ensemble/),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="steps",&
            description="The number of MD steps to perform",&
            usage="steps 100", default_i_val=3,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,name="timestep",&
            description="The length of an integration step (in case RESPA the large TIMESTEP)",&
            usage="timestep 1.0",default_r_val=cp_unit_to_cp2k(value=0.5_dp,unit_str="fs",error=error),&
            unit_str="fs",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,name="step_start_val",&
            description="The starting step value for the MD",usage="step_start_val <integer>",&
            default_i_val=0,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,name="time_start_val",&
            description="The starting timer value for the MD",&
            usage="time_start_val <real>",default_r_val=cp_unit_to_cp2k(value=0.0_dp,unit_str="fs",error=error),&
            unit_str="fs",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,name="econs_start_val",&
            description="The starting  value of the conserved quantity",&
            usage="econs_start_val <real>",default_r_val=0.0_dp,&
            unit_str="hartree",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,name="temperature",&
            description="The temperature in K used to initialize "//&
            "the velocities with init and pos restart, and in the NPT/NVT simulations",&
            usage="temperature 325.0",default_r_val=cp_unit_to_cp2k(value=300.0_dp,unit_str="K",error=error),&
            unit_str="K",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="temp_tol",&
            variants=s2a("temp_to","temperature_tolerance"),&
            description="The maximum accepted deviation of the (global) temperature"//&
            "from the desired target temperature before a rescaling of the velocites "//&
            "is performed. If it is 0 no rescaling is performed. NOTE: This keyword is "//&
            "obsolescent; Using a CSVR thermostat with a short timeconstant is "//&
            "recommended as a better alternative.", &
            usage="temp_tol 0.0", default_r_val=0.0_dp, unit_str='K', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="temp_kind",&
            description="Compute the temperature per each kind separately",&
            usage="temp_kind LOGICAL",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="scale_temp_kind",&
            description="When necessary rescale the temperature per each kind separately",&
            usage="scale_temp_kind LOGICAL",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="comvel_tol",&
            description="The maximum accepted velocity of the center of mass. "//&
            "With Shell-Model, comvel may drift if MD%THERMOSTAT%REGION /= GLOBAL ", &
            usage="comvel_tol 0.1", type_of_var=real_t, n_var=1, unit_str="bohr*au_t^-1",&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="angvel_tol",&
            description="The maximum accepted angular velocity. This option is ignored "//&
            "when the system is periodic. Removes the components of the velocities that"//&
            "project on the external rotational degrees of freedom.",&
            usage="angvel_tol 0.1", type_of_var=real_t, n_var=1, unit_str="bohr*au_t^-1",&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="angvel_zero",&
            description="Set the initial angular velocity to zero. This option is ignored "//&
            "when the system is periodic or when initial velocities are defined. Technically, "//&
            "the part of the random initial velocities that projects on the external "//&
            "rotational degrees of freedom is subtracted.",&
            usage="angvel_zero LOGICAL",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ANNEALING",&
            description="Specifies the rescaling factor for annealing velocities. "//&
            "Automatically enables the annealing procedure. This scheme works only for ensembles "//&
            "that do not have thermostats on particles.", &
            usage="annealing <REAL>", default_r_val=1.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ANNEALING_CELL",&
            description="Specifies the rescaling factor for annealing velocities of the CELL "//&
            "Automatically enables the annealing procedure for the CELL. This scheme works only "//&
            "for ensambles that do not have thermostat on CELLS velocities.", &
            usage="ANNEALING_CELL <REAL>", default_r_val=1.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DISPLACEMENT_TOL",&
            description="This keyword sets a maximum atomic displacement "//&
            " in each Cartesian direction."//&
            "The maximum velocity is evaluated and if it is too large to remain"//&
            "within the assigned limit, the time step is rescaled accordingly,"//&
            "and the first half step of the velocity verlet is repeated.",&
            usage="DISPLACEMENT_TOL <REAL>", default_r_val=100.0_dp,&
            unit_str='angstrom',  error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_langevin_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_msst_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_barostat_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_thermostat_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_respa_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_shell_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_adiabatic_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_softening_section(subsection, error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_reftraj_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_avgs_section(subsection, error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection, error=error)

       CALL create_thermal_region_section(subsection, error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection, error=error)

       CALL create_md_print_section(subsection, error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection, error=error)

       CALL create_cascade_section(subsection, error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection, error=error)
    END IF

  END SUBROUTINE create_md_section

! *****************************************************************************
!> \brief Defines AVERAGES section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_avgs_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_avgs_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="Averages",&
            description="Controls the calculation of the averages during an MD run.",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)
       NULLIFY(keyword, print_key, subsection)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="Controls the calculations of the averages.",&
            usage="&AVERAGES T",default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ACQUISITION_START_TIME",&
            description="Setup up the simulation time when the acquisition process to compute "//&
            " averages is started.",&
            usage="ACQUISITION_START_TIME <REAL>",&
            default_r_val=0.0_dp, unit_str='fs',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVERAGE_COLVAR",&
            description="Switch for computing the averages of COLVARs.", &
            usage="AVERAGE_COLVAR <LOGICAL>", default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL cp_print_key_section_create(print_key,"PRINT_AVERAGES",&
            description="Controls the output the averaged quantities",&
            print_level=debug_print_level+1, common_iter_levels=1,&
            filename="",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL create_avgs_restart_section(subsection, error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection, error=error)
    END IF
  END SUBROUTINE create_avgs_section

! *****************************************************************************
!> \brief Defines the AVERAGES RESTART section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_avgs_restart_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_avgs_restart_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="RESTART_AVERAGES",&
            description="Stores information for restarting averages.",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ITIMES_START",&
            description="TIME STEP starting the evaluation of averages",&
            usage="ITIMES_START <INTEGER>",type_of_var=integer_t,n_var=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVECPU",&
            description="CPU average",usage="AVECPU <REAL>",&
            type_of_var=real_t,n_var=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVEHUGONIOT",&
            description="HUGONIOT average",usage="AVEHUGONIOT <REAL>",&
            type_of_var=real_t,n_var=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVETEMP_BARO",&
            description="BAROSTAT TEMPERATURE average",usage="AVETEMP_BARO <REAL>",&
            type_of_var=real_t,n_var=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVEPOT",&
            description="POTENTIAL ENERGY average",usage="AVEPOT <REAL>",&
            type_of_var=real_t,n_var=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVEKIN",&
            description="KINETIC ENERGY average",usage="AVEKIN <REAL>",&
            type_of_var=real_t,n_var=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVETEMP",&
            description="TEMPERATURE average",usage="AVETEMP <REAL>",&
            type_of_var=real_t,n_var=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVEKIN_QM",&
            description="QM KINETIC ENERGY average in QMMM runs",usage="AVEKIN_QM <REAL>",&
            type_of_var=real_t,n_var=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVETEMP_QM",&
            description="QM TEMPERATURE average in QMMM runs",usage="AVETEMP_QM <REAL>",&
            type_of_var=real_t,n_var=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVEVOL",&
            description="VOLUME average",usage="AVEVOL <REAL>",&
            type_of_var=real_t,n_var=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVECELL_A",&
            description="CELL VECTOR A average",usage="AVECELL_A <REAL>",&
            type_of_var=real_t,n_var=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVECELL_B",&
            description="CELL VECTOR B average",usage="AVECELL_B <REAL>",&
            type_of_var=real_t,n_var=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVECELL_C",&
            description="CELL VECTOR C average",usage="AVECELL_C <REAL>",&
            type_of_var=real_t,n_var=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVEALPHA",&
            description="ALPHA cell angle average",usage="AVEALPHA <REAL>",&
            type_of_var=real_t,n_var=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVEBETA",&
            description="BETA cell angle average",usage="AVEBETA <REAL>",&
            type_of_var=real_t,n_var=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVEGAMMA",&
            description="GAMMA cell angle average",usage="AVEGAMMA <REAL>",&
            type_of_var=real_t,n_var=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVE_ECONS",&
            description="CONSTANT ENERGY average",usage="AVE_ECONS <REAL>",&
            type_of_var=real_t,n_var=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVE_PRESS",&
            description="PRESSURE average",usage="AVE_PRESS <REAL>",&
            type_of_var=real_t,n_var=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVE_PXX",&
            description="P_{XX} average",usage="AVE_PXX <REAL>",&
            type_of_var=real_t,n_var=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVE_PV_VIR",&
            description="PV VIRIAL average",usage="AVE_PV_VIR <REAL> .. <REAL>",&
            type_of_var=real_t, n_var=9, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVE_PV_TOT",&
            description="PV TOTAL average",usage="AVE_PV_TOT <REAL> .. <REAL>",&
            type_of_var=real_t, n_var=9, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVE_PV_KIN",&
            description="PV KINETIC average",usage="AVE_PV_KIN <REAL> .. <REAL>",&
            type_of_var=real_t, n_var=9, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVE_PV_CNSTR",&
            description="PV CONSTRAINTS average",usage="AVE_PV_CNSTR <REAL> .. <REAL>",&
            type_of_var=real_t, n_var=9, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVE_PV_XC",&
            description="PV XC average",usage="AVE_PV_XC <REAL> .. <REAL>",&
            type_of_var=real_t, n_var=9, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVE_PV_FOCK_4C",&
            description="PV XC average",usage="AVE_PV_FOCK_4C <REAL> .. <REAL>",&
            type_of_var=real_t, n_var=9, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVE_COLVARS",&
            description="COLVARS averages",usage="AVE_COLVARS <REAL> .. <REAL>",&
            type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AVE_MMATRIX",&
            description="METRIC TENSOR averages",usage="AVE_MMATRIX <REAL> .. <REAL>",&
            type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF
  END SUBROUTINE create_avgs_restart_section

! *****************************************************************************
!> \brief Defines LANGEVIN section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_langevin_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_langevin_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="Langevin",&
            description="Controls the set of parameters to run a Langevin MD",&
            citations=(/Ricci2003/),&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="gamma",&
            description="Gamma parameter for the Langevin dynamics (LD)",&
            usage="gamma 0.001",&
            default_r_val=0.0_dp, unit_str='fs^-1',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="Noisy_Gamma",&
            variants=(/"NoisyGamma"/), &
            description="Imaginary Langevin Friction term for LD with noisy forces.",&
            citations=(/Kuhne2007/),&
            usage="Noisy_Gamma 4.0E-5", default_r_val=0.0_dp, unit_str='fs^-1', error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="Shadow_Gamma",&
            variants=(/"ShadowGamma"/), &
            description="Shadow Langevin Friction term for LD with noisy forces in order to adjust Noisy_Gamma.", &
            citations=(/Kuhne2007/),&
            usage="Shadow_Gamma 0.001", default_r_val=0.0_dp, unit_str='fs^-1', error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)
    END IF
  END SUBROUTINE create_langevin_section

! *****************************************************************************
!> \brief Defines print section for MD
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_md_print_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_md_print_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="print",&
            description="Controls the printing properties during an MD run",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)
       NULLIFY(print_key, keyword)

       CALL keyword_create(keyword, name="FORCE_LAST",&
            description="Print the output and restart file if walltime is reached or "//&
           "if an external EXIT command is given. It still requires the keyword LAST "//&
           "to be present for the specific print key (in case the last step should not "//&
           "match with the print_key iteration number).",&
            usage="FORCE_LAST LOGICAL",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL cp_print_key_section_create(print_key,"ENERGY",&
            description="Controls the output the ener file",&
            print_level=low_print_level, common_iter_levels=1,&
            filename="",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SHELL_ENERGY",&
            description="Controls the output of the shell-energy file (only if shell-model)",&
            print_level=medium_print_level, common_iter_levels=1,&
            filename="",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"TEMP_KIND",&
            description="Controls the output of the temperature"//&
            " computed separately for each kind",&
            print_level=high_print_level, common_iter_levels=1,&
            filename="",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"TEMP_SHELL_KIND",&
            description="Controls the output of the temperature of the"//&
            " shell-core motion computed separately for each kind",&
            print_level=high_print_level, common_iter_levels=1,&
            filename="",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"CENTER_OF_MASS",&
            description="Controls the printing of COM velocity during an MD", &
            print_level=medium_print_level,common_iter_levels=1,&
            filename="__STD_OUT__",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"COEFFICIENTS",&
            description="Controls the printing of coefficients during an MD run.", &
            print_level=medium_print_level,common_iter_levels=1,&
            filename="",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"ROTATIONAL_INFO",&
            description="Controls the printing basic info during the calculation of the "//&
            "translational/rotational degrees of freedom.", print_level=low_print_level,&
            add_last=add_last_numeric,filename="__STD_OUT__",error=error)
       CALL keyword_create(keyword, name="COORDINATES",&
            description="Prints atomic coordinates in the standard orientation. "//&
            "Coordinates are not affected during the calculation.",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="Controls the printing of basic and summary information during the"//&
            " Molecular Dynamics", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)
    END IF
  END SUBROUTINE create_md_print_section

! *****************************************************************************
!> \brief Defines parameters for RESPA integration scheme
!> \param section will contain the coeff section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_respa_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_respa_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="RESPA",&
            description="Multiple timestep integration based on RESPA (implemented for NVE only)."//&
            " RESPA exploits multiple force_eval."//&
            " In this case the order of the force_eval maps "//&
            " the order of the respa shells from the slowest to the fastest force evaluation."//&
            " If force_evals share the same subsys, it's enough then to specify the "//&
            " subsys in the force_eval corresponding at the first index in the multiple_force_eval list."//&
            " Can be used to speedup classical and ab initio MD simulations.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            citations=(/Tuckerman1992,Guidon2008/),&
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword, name="FREQUENCY",&
            description="The number of reference MD steps between two RESPA corrections.",&
            usage="FREQUENCY <INTEGER>", default_i_val=5,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_respa_section

! *****************************************************************************
!> \brief Defines parameters for REFTRAJ analysis
!> \param section will contain the coeff section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_reftraj_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_reftraj_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="REFTRAJ",&
            description="Loads an external trajectory file and performs analysis on the"//&
            " loaded snapshots.",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword, print_key, subsection)

       CALL keyword_create(keyword, name="TRAJ_FILE_NAME",&
            description="Specify the filename where the trajectory is stored.",repeats=.FALSE.,&
            usage="TRAJ_FILE_NAME <CHARACTER>", default_lc_val="reftraj.xyz", &
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CELL_FILE_NAME",&
            description="Specify the filename where the cell is stored "//&
            "(for trajectories generated within variable cell ensembles).",repeats=.FALSE.,&
            usage="CELL_FILE_NAME <CHARACTER>", default_lc_val="reftraj.cell", &
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="VARIABLE_VOLUME",&
            description="Enables the possibility to read a CELL file with information on the CELL size during the MD.",&
            repeats=.FALSE., default_l_val=.FALSE.,lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FIRST_SNAPSHOT",&
            description="Index of the snapshot stored in the trajectory file "//&
            "from which to start a REFTRAJ run",&
            repeats=.FALSE., usage="FIRST_SNAPSHOT <INTEGER>", default_i_val=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LAST_SNAPSHOT",&
            description="Index of the last snapshot stored in the trajectory file that is read along a REFTRAJ run",&
            repeats=.FALSE., usage="LAST_SNAPSHOT", default_i_val=0, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="STRIDE",&
            description=" Stride in number of snapshot for the  reftraj analysis",&
            repeats=.FALSE., usage="STRIDE", default_i_val=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="eval_energy_forces",&
            description="Evaluate energy and forces for each retrieved snapshot during a REFTRAJ run",&
            repeats=.FALSE., default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_msd_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="print",&
       description="The section that controls the output of a reftraj run",&
       n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
       error=error)

       NULLIFY(print_key)
         CALL cp_print_key_section_create(print_key,"msd_kind",&
              description="Controls the output of msd per kind", &
              print_level=low_print_level, common_iter_levels=1,&
              filename="",error=error)
         CALL section_add_subsection(subsection,print_key,error=error)
         CALL section_release(print_key,error=error)

         CALL cp_print_key_section_create(print_key,"msd_molecule",&
              description="Controls the output of msd per molecule kind", &
              print_level=low_print_level, common_iter_levels=1,&
              filename="",error=error)
         CALL section_add_subsection(subsection,print_key,error=error)
         CALL section_release(print_key,error=error)

         CALL cp_print_key_section_create(print_key,"displaced_atom",&
              description="Controls the output of index and dislacement of "//&
              "atoms that moved away from the initial position of more than a"//&
              "given distance (see msd%disp_tol)", &
              print_level=low_print_level, common_iter_levels=1,&
              filename="",error=error)
         CALL section_add_subsection(subsection,print_key,error=error)
         CALL section_release(print_key,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection, error=error)

    END IF
  END SUBROUTINE create_reftraj_section

! *****************************************************************************
!> \brief Defines parameters for MSD calculation along a REFTRAJ analysis
!> \param section will contain the coeff section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author MI
! *****************************************************************************
  SUBROUTINE create_msd_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_msd_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="MSD",&
            description="Loads an external trajectory file and performs analysis on the"//&
            " loaded snapshots.",&
            n_keywords=3, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="controls the activation of core-level spectroscopy simulations",&
            usage="&MSD T",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="REF0_FILENAME",&
            description="Specify the filename where the initial reference configuration is stored.",&
            repeats=.FALSE., usage="REF0_FILENAME <CHARACTER>", default_lc_val="", &
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="MSD_PER_KIND",&
            description="Set up the calculation of the MSD for each atomic kind",&
            usage="MSD_PER_KIND <LOGICAL>",repeats=.FALSE.,&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="MSD_PER_MOLKIND",&
            description="Set up the calculation of the MSD for each molecule kind."//&
            "The position of the center of mass of the molecule is considered.",&
            usage="MSD_PER_MOLKIND <LOGICAL>",repeats=.FALSE.,&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="MSD_PER_REGION",&
            description="Set up the calculation of the MSD for each defined region.",&
            usage="MSD_PER_REGION <LOGICAL>",repeats=.FALSE.,&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_region_section(subsection,"MSD calculation",error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL keyword_create(keyword,name="DISPLACED_ATOM",&
            description="Identify the atoms that moved from their initial"//&
            "position of a distance larger than a given tolerance (see msd%displacement_tol).",&
            usage="DISPLACED_ATOM <LOGICAL>",repeats=.FALSE.,&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="displacement_tol",&
              description="Lower limit to define displaced atoms",&
              usage="DISPLACEMENT_TOL real",&
              default_r_val=0._dp, n_var=1, unit_str='bohr', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_msd_section

! *****************************************************************************
!> \param section will contain the coeff section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_msst_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_msst_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="msst",&
            description="Parameters for Multi-Scale Shock Technique (MSST) "//&
            "which simulate the effect of a steady planar shock on a unit cell. "//&
            "Reed et. al. Physical Review Letters 90, 235503 (2003).",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword )
       CALL keyword_create(keyword, name="PRESSURE",&
            description="Initial pressure",&
            usage="PRESSURE real",&
            default_r_val=0._dp, n_var=1, unit_str='bar', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ENERGY",&
            description="Initial energy",&
            usage="ENERGY real",&
            default_r_val=0._dp, n_var=1, unit_str='hartree',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="VOLUME",&
            description="Initial volume",&
            usage="VOLUME real",&
            default_r_val=0._dp, n_var=1, unit_str='angstrom^3',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CMASS",&
            description="Effective cell mass",&
            usage="CMASS real",&
            default_r_val=0._dp, n_var=1,unit_str='au_m',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="VSHOCK",variants=(/"V_SHOCK"/),&
            description="Velocity shock",&
            usage="VSHOCK real",&
            default_r_val=0._dp, n_var=1,unit_str='m/s',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="GAMMA",&
            description="Damping coefficient for cell volume",&
            usage="GAMMA real",&
            unit_str='fs^-1', &
            default_r_val=0.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_msst_section

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
         repeats=.TRUE., usage="<REAL> <REAL> .. <REAL>",&
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
         "helium atom",repeats=.TRUE.,usage="<INT> <INT> .. <INT>",&
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
         usage="<REAL> <REAL> .. <REAL>", type_of_var=real_t, &
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

! *****************************************************************************
!> \brief creates the FORMAT keyword
!> \param section will contain the pint section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino 10.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE add_format_keyword(keyword, section, pos, description, error)
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section
    LOGICAL, INTENT(IN)                      :: pos
    CHARACTER(LEN=*), INTENT(IN)             :: description
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'add_format_keyword', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(keyword),cp_failure_level,routineP,error,failure)

    IF (pos) THEN

       CALL keyword_create(keyword, name="FORMAT",&
            description=description,usage="FORMAT (ATOMIC|DCD|PDB|XMOL|XYZ)",&
            default_i_val=dump_xmol,&
            enum_c_vals=s2a("ATOMIC","DCD","PDB","XMOL","XYZ"),&
            enum_i_vals= (/dump_atomic,dump_dcd,dump_pdb,dump_xmol,dump_xmol/),&
            enum_desc=s2a("Write only the coordinates X,Y,Z without element symbols to a formatted file",&
                          "Write the coordinates (no element labels) and the cell information to a binary file",&
                          "Write the atomic information in PDB format to a formatted file",&
                          "Mostly known as XYZ format, provides in a formatted file: element_symbol X Y Z",&
                          "Alias name for XMOL"),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CHARGE_OCCUP",&
            variants=(/"CHARGE_O"/),&
            description="Write the MM charges to the OCCUP field of the PDB file",&
            usage="CHARGE_OCCUP logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CHARGE_BETA",&
            variants=(/"CHARGE_B"/),&
            description="Write the MM charges to the BETA field of the PDB file",&
            usage="CHARGE_BETA logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CHARGE_EXTENDED",&
            description="Write the MM charges to the very last field of the PDB file (starting from column 81)",&
            usage="CHARGE_EXTENDED logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    ELSE

       CALL keyword_create(keyword, name="FORMAT",&
            description=description,usage="FORMAT (ATOMIC|DCD|XMOL|XYZ)",&
            default_i_val=dump_xmol,&
            enum_c_vals=s2a("ATOMIC","DCD","XMOL","XYZ"),&
            enum_i_vals= (/dump_atomic,dump_dcd,dump_xmol,dump_xmol/),&
            enum_desc=s2a("Write only the coordinates X,Y,Z without element symbols to a formatted file",&
                          "Write the coordinates (no element labels) and the cell information to a binary file",&
                          "Mostly known as XYZ format, provides in a formatted file: element_symbol X Y Z",&
                          "Alias name for XMOL"),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE add_format_keyword

! *****************************************************************************
!> \brief creates the section for a BAND run
!> \param section will contain the pint section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino 09.2006 [tlaino]
! *****************************************************************************
  SUBROUTINE create_band_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_band_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection, &
                                                subsubsection

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="band",&
            description="The section that controls a BAND run",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            citations=(/Elber1987,Jonsson1998,Jonsson2000_1,Jonsson2000_2,Wales2004/),&
            error=error)
       NULLIFY(keyword, print_key, subsection, subsubsection)

       CALL keyword_create(keyword, name="NPROC_REP",&
            description="Specify the number of processors to be used per replica "//&
            "environment (for parallel runs)",&
            default_i_val=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PROC_DIST_TYPE",&
            description="Specify the topology of the mapping of processors into replicas.",&
            usage="PROC_DIST_TYPE (INTERLEAVED|BLOCKED)",&
            enum_c_vals=s2a("INTERLEAVED",&
                            "BLOCKED"),&
            enum_desc=s2a( "Interleaved distribution",&
                           "Blocked distribution"),&
            enum_i_vals=(/do_rep_interleaved,do_rep_blocked/),&
            default_i_val=do_rep_blocked, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

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
            enum_desc=s2a( "Bisection NEB",&
                           "Improved tangent NEB",&
                           "Climbing Image NEB",&
                           "Doubly NEB",&
                           "String Method",&
                           "Elastic Band (Hamiltonian formulation)"),&
            enum_i_vals=(/do_b_neb,do_it_neb,do_ci_neb,do_d_neb,do_sm,do_eb/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NUMBER_OF_REPLICA",&
            description="Specify the number of Replica to use in the BAND",&
            default_i_val=10, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="USE_COLVARS",&
            description="Uses a version of the band scheme projected in a subspace of colvars.",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

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
            enum_i_vals=(/pot_neb_full,pot_neb_fe,pot_neb_me/),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ROTATE_FRAMES",&
            description="Compute at each BAND step the RMSD and rotate the frames in order"//&
            " to minimize it.",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ALIGN_FRAMES",&
            description="Enables the alignment of the frames at the beginning of a BAND calculation. "//&
            "This keyword does not affect the rotation of the replicas during a BAND calculation.",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="K_SPRING",&
            variants=(/"K"/),&
            description="Specify the value of the spring constant",&
            default_r_val=0.02_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Convergence_control
       CALL section_create(subsection,name="CONVERGENCE_CONTROL",&
            description="Setup parameters to control the convergence criteria for BAND",&
            repeats=.FALSE., required=.FALSE., error=error)
       CALL keyword_create(keyword, name="MAX_DR",&
            description="Tolerance on the maximum value of the displacement on the BAND.",&
            usage="MAX_DR {real}",&
            default_r_val=0.0002_dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_FORCE",&
            description="Tolerance on the maximum value of Forces on the BAND.",&
            usage="MAX_FORCE {real}",&
            default_r_val=0.00045_dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMS_DR",&
            description="Tolerance on RMS displacements on the BAND.",&
            usage="RMS_DR {real}",&
            default_r_val=0.0001_dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMS_FORCE",&
            description="Tolerance on RMS Forces on the BAND.",&
            usage="RMS_FORCE {real}",&
            default_r_val=0.00030_dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection, subsubsection)
       ! CI-NEB section
       CALL section_create(subsection,name="CI_NEB",&
            description="Controls parameters for CI-NEB type calculation only.",&
            repeats=.FALSE., required=.FALSE., error=error)
       CALL keyword_create(keyword, name="NSTEPS_IT",&
            description="Specify the number of steps of IT-NEB to perform before "//&
            "switching on the CI algorithm",&
            default_i_val=5, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       ! String Method section
       CALL section_create(subsection,name="STRING_METHOD",&
            description="Controls parameters for String Method type calculation only.",&
            repeats=.FALSE., required=.FALSE., error=error)

       CALL keyword_create(keyword, name="SPLINE_ORDER",&
            description="Specify the oder of the spline used in the String Method.",&
            default_i_val=1, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="SMOOTHING",&
            description="Smoothing parameter for the reparametrization of the frames.",&
            default_r_val=0.2_dp, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       ! Optimization section
       CALL section_create(subsection,name="optimize_band",&
            description="Specify the optimization method for the band",&
            repeats=.TRUE., required=.TRUE., error=error)
       CALL create_opt_band_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       ! replica section: to specify coordinates and velocities (possibly) of the
       ! different replica used in the BAND
       CALL section_create(subsection,name="replica",&
            description="Specify coordinates and velocities (possibly) of the replica",&
            repeats=.TRUE., required=.TRUE., error=error)
       ! Colvar
       CALL keyword_create(keyword, name="COLLECTIVE",&
            description="Specifies the value of the collective variables used in the projected"//&
            " BAND method. The order of the values is the order of the COLLECTIVE section in the"//&
            " constraints/restraints section",&
            usage="COLLECTIVE {real} .. {real}",&
            type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       ! Coordinates read through an external file
       CALL keyword_create(keyword, name="COORD_FILE_NAME",&
            description="Name of the xyz file with coordinates (alternative to &COORD section)",&
            usage="COORD_FILE_NAME <CHAR>",&
            default_lc_val="",error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       ! Coordinates and velocities
       CALL create_coord_section(subsubsection,"BAND",error=error)
       CALL section_add_subsection(subsection,subsubsection,error=error)
       CALL section_release(subsubsection,error=error)
       CALL create_velocity_section(subsubsection,"BAND",error=error)
       CALL section_add_subsection(subsection,subsubsection,error=error)
       CALL section_release(subsubsection,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       ! Print key section
       CALL cp_print_key_section_create(print_key,"program_run_info",&
            description="Controls the printing basic info about the BAND run", &
            print_level=medium_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)

       CALL keyword_create(keyword, name="INITIAL_CONFIGURATION_INFO",&
            description="Print information for the setup of the initial configuration.",&
            usage="INITIAL_CONFIGURATION_INFO <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"convergence_info",&
            description="Controls the printing of the convergence criteria during a BAND run", &
            print_level=medium_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"replica_info",&
            description="Controls the printing of each replica info during a BAND run", &
            print_level=medium_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"ENERGY",&
            description="Controls the printing of the ENER file in a BAND run", &
            print_level=low_print_level, common_iter_levels=1,&
            filename="",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"BANNER",&
            description="Controls the printing of the BAND banner", &
            print_level=low_print_level, common_iter_levels=1,&
            filename="__STD_OUT__",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)
    END IF
  END SUBROUTINE create_band_section

! *****************************************************************************
!> \brief creates the optimization section for a BAND run
!> \param section will contain the pint section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino 02.2007 [tlaino]
! *****************************************************************************
  SUBROUTINE create_opt_band_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_opt_band_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection, &
                                                subsubsection

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(keyword, print_key, subsection, subsubsection)

       CALL keyword_create(keyword, name="OPT_TYPE",&
            description="Specifies the type optimizer used for the band",&
            usage="OPT_TYPE (MD|DIIS)",&
            default_i_val=band_diis_opt,&
            enum_c_vals=s2a("MD",&
                            "DIIS"),&
            enum_desc=s2a( "Molecular Dynamics based Optimizer",&
                           "Coupled Steepest Descent / DIIS "),&
            enum_i_vals=(/band_md_opt,band_diis_opt/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="OPTIMIZE_END_POINTS",&
            description="Performs also an optimization of the end points of the band.",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! MD optimization section
       CALL section_create(subsection,name="MD",&
            description="Activate the MD based optimization procedure for BAND",&
            repeats=.FALSE., required=.TRUE., error=error)

       CALL keyword_create(keyword, name="MAX_STEPS",&
            description="Specify the maximum number of MD steps",&
            default_i_val=100, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(&
            keyword=keyword,&
            name="timestep",&
            description="The length of an integration step",&
            usage="timestep 1.0",&
            default_r_val=cp_unit_to_cp2k(value=0.5_dp,&
            unit_str="fs",&
            error=error),&
            unit_str="fs",&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TEMPERATURE",&
            description="Specify the initial temperature",&
            default_r_val=cp_unit_to_cp2k(value=0.0_dp,&
            unit_str="K",&
            error=error),&
            unit_str="K",&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Temp_control
       CALL section_create(subsubsection,name="TEMP_CONTROL",&
            description="Setup parameters to control the temperature during a BAND MD run.",&
            repeats=.FALSE., required=.TRUE., error=error)
       CALL keyword_create(keyword, name="TEMPERATURE",&
            description="Specify the target temperature",&
            type_of_var=real_t,unit_str="K",error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TEMP_TOL",&
            description="Specify the tolerance on the temperature for rescaling",&
            default_r_val=cp_unit_to_cp2k(value=0.0_dp,&
            unit_str="K",&
            error=error),&
            unit_str="K",&
            error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TEMP_TOL_STEPS",&
            description="Specify the number of steps to apply a temperature control",&
            default_i_val=0, error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(subsection, subsubsection, error=error)
       CALL section_release(subsubsection,error=error)

       ! Vel_control
       CALL section_create(subsubsection,name="VEL_CONTROL",&
            description="Setup parameters to control the velocity during a BAND MD run.",&
            repeats=.FALSE., required=.TRUE., error=error)
       CALL keyword_create(keyword, name="ANNEALING",&
            description="Specify the annealing coefficient",&
            default_r_val=1.0_dp, error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="PROJ_VELOCITY_VERLET",&
            description="Uses a Projected Velocity Verlet instead of a normal Velocity Verlet."//&
            " Every time the cosine between velocities and forces is < 0 velocities are"//&
            " zeroed.",&
            usage="PROJ_VELOCITY_VERLET <LOGICAL>",&
            default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="SD_LIKE",&
            description="Zeros velocity at each MD step emulating a steepest descent like"//&
            "(SD_LIKE) approach",&
            usage="SD_LIKE <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(subsection, subsubsection, error=error)
       CALL section_release(subsubsection,error=error)
       ! End of MD
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       ! DIIS optimization section
       CALL section_create(subsection,name="DIIS",&
            description="Activate the DIIS based optimization procedure for BAND",&
            repeats=.FALSE., required=.TRUE., error=error)

       CALL keyword_create(keyword, name="MAX_SD_STEPS",&
            description="Specify the maximum number of SD steps to perform"//&
            " before switching on DIIS (the minimum number will always be equal to N_DIIS).",&
            default_i_val=1, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_STEPS",&
            description="Specify the maximum number of optimization steps",&
            default_i_val=100, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="N_DIIS",&
            variants=(/"NDIIS"/),&
            description="Number of history vectors to be used with DIIS",&
            usage="N_DIIS 4",&
            default_i_val=7,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="STEPSIZE",&
            description="Initial stepsize used for the line search, sometimes this parameter"//&
            "can be reduced to stablize DIIS",&
            usage="STEPSIZE <REAL>",&
            default_r_val=1.0_dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_STEPSIZE",&
            description="Maximum stepsize used for the line search, sometimes this parameter"//&
            "can be reduced to stablize the LS for particularly difficult initial geometries",&
            usage="MAX_STEPSIZE <REAL>",&
            default_r_val=2.0_dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NP_LS",&
            description="Number of points used in the line search SD.",&
            usage="NP_LS <INTEGER>",&
            default_i_val=2,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NO_LS",&
            description="Does not perform LS during SD. Useful in combination with a proper STEPSIZE"//&
            " for particularly out of equilibrium starting geometries.",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CHECK_DIIS",&
            description="Performes a series of checks on the DIIS solution in order to accept the DIIS step."//&
            " If set to .FALSE. the only check performed is that the angle between the DIIS solution and the"//&
            " reference vector is less than Pi/2. Can be useful if many DIIS steps are rejected.",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL cp_print_key_section_create(print_key,"diis_info",&
            description="Controls the printing of diis info during a BAND run", &
            print_level=high_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
    END IF
  END SUBROUTINE create_opt_band_section

! *****************************************************************************
!> \brief section will contain some parameters for the adiabatic dynamics
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE create_adiabatic_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_adiabatic_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: thermo_fast_section, &
                                                thermo_slow_section

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="ADIABATIC_DYNAMICS",&
            description="Parameters used in canonical adiabatic free energy sampling (CAFES).",&
            n_keywords=5, n_subsections=2, repeats=.FALSE., required=.FALSE.,&
            citations=(/VandeVondele2002/),&
            error=error)

       NULLIFY(keyword,thermo_fast_section,thermo_slow_section)

       CALL keyword_create(keyword=keyword,name="temp_fast",&
            description="Temperature in K used to control "//&
            "the fast degrees of freedom ",&
            usage="temp_fast 5.0",&
            default_r_val=cp_unit_to_cp2k(value=0.0_dp,unit_str="K",error=error),&
            unit_str="K", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,name="temp_slow",&
            description="Temperature in K used to control "//&
            "the slow degrees of freedom ",&
            usage="temp_slow 5.0",&
            default_r_val=cp_unit_to_cp2k(value=0.0_dp,unit_str="K",error=error),&
            unit_str="K", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="temp_tol_fast",&
            description="Maximum accepted temperature deviation"//&
            " from the expected value, for the fast motion."//&
            "If 0, no rescaling is performed", &
            usage="temp_tol 0.0", default_r_val=0.0_dp, unit_str='K', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="temp_tol_slow",&
            description="Maximum accepted temperature deviation"//&
            " from the expected value, for the slow motion."//&
            "If 0, no rescaling is performed", &
            usage="temp_tol 0.0", default_r_val=0.0_dp, unit_str='K', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="n_resp_fast",&
            description="number of respa steps for fast degrees of freedom",&
            repeats=.FALSE., default_i_val=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL create_thermo_fast_section(thermo_fast_section, error=error)
       CALL section_add_subsection(section, thermo_fast_section, error=error)
       CALL section_release(thermo_fast_section,error=error)

       CALL create_thermo_slow_section(thermo_slow_section, error=error)
       CALL section_add_subsection(section, thermo_slow_section, error=error)
       CALL section_release(thermo_slow_section,error=error)

    END IF
  END SUBROUTINE create_adiabatic_section


! *****************************************************************************
!> \brief section will contain parameters for the velocity softening
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE create_softening_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_softening_section', &
      routineP = moduleN//':'//routineN

    TYPE(keyword_type), POINTER              :: keyword

    CALL section_create(section,name="VELOCITY_SOFTENING",&
         description="A method to initialize the velocities along low-curvature "&
         //"directions in order to favors MD trajectories to cross rapidly over "&
         //"small energy barriers into neighboring basins. "&
         //"In each iteration the forces are calculated at a point y, which "&
         //"is slightly displaced from the current positions x in the direction "&
         //"of the original velocities v. The velocities are then updated with "&
         //"the force component F_t, which is perpendicular to N. "&
         //"N = v / |v|;  y = x + delta * N;  F_t = F(y) - <F(y) | N> * N;  "&
         //"v' = v + alpha * F_t", error=error)

    NULLIFY(keyword)
    CALL keyword_create(keyword=keyword, name="STEPS",&
         description="Number of softening iterations performed. "&
         //"Typical values are around 40 steps.",&
         default_i_val=0, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword=keyword, name="DELTA",&
         description="Displacement used to obtain y.",&
         default_r_val=0.1_dp, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword=keyword, name="ALPHA",&
         description="Mixing factor used for updating velocities.",&
         default_r_val=0.15_dp, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

  END SUBROUTINE create_softening_section


! *****************************************************************************
!> \brief section will contain some parameters for the shells dynamics
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE create_shell_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_shell_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: thermo_section

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="shell",&
            description="Parameters of shell model in adiabatic dynamics.",&
            n_keywords=4, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword,thermo_section)

       CALL keyword_create(keyword=keyword,name="temperature",&
            description="Temperature in K used to control "//&
            "the internal velocities of the core-shell motion ",&
            usage="temperature 5.0",&
            default_r_val=cp_unit_to_cp2k(value=0.0_dp,unit_str="K",error=error),&
            unit_str="K", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="temp_tol",&
            description="Maximum accepted temperature deviation"//&
            " from the expected value, for the internal core-shell motion."//&
            "If 0, no rescaling is performed", &
            usage="temp_tol 0.0", default_r_val=0.0_dp, unit_str='K', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="nose_particle",&
            description="If nvt or npt, the core and shell velocities are controlled "//&
            "by the same thermostat used for the particle. This might favour heat exchange "//&
            "and additional rescaling of the internal core-shell velocity is needed (TEMP_TOL)",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DISPLACEMENT_SHELL_TOL",&
            description="This keyword sets a maximum variation of the shell "//&
            "core distance in each Cartesian direction."//&
            "The maximum internal core-shell velocity is evaluated and"//&
            " if it is too large to remain"//&
            "within the assigned limit, the time step is rescaled accordingly,"//&
            "and the first half step of the velocity verlet is repeated.",&
            usage="DISPLACEMENT_SHELL_TOL <REAL>", default_r_val=100.0_dp,&
            unit_str='angstrom',  error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_thermostat_section(thermo_section, error=error)
       CALL section_add_subsection(section, thermo_section, error=error)
       CALL section_release(thermo_section,error=error)

    END IF
  END SUBROUTINE create_shell_section

! *****************************************************************************
!> \brief input section used to define regions with different temperature
!>        initialization and control
! *****************************************************************************
  SUBROUTINE create_thermal_region_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_thermal_region_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, region_section, &
                                                subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="thermal_region",&
            description="Define regions where different initialization and control "//&
            "of the temperature is used.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)


       NULLIFY(region_section)
       NULLIFY(keyword,subsection)

       CALL keyword_create(keyword=keyword,name="force_rescaling",&
            description="Control the rescaling ot the velocities in all the regions, "//&
            " according to the temperature assigned to each reagion, when RESTART_VELOCITY in EXT_RESTART is active.",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_create(region_section,name="DEFINE_REGION",&
            description="This section provides the possibility to define arbitrary region ",&
            n_keywords=3, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword, name="LIST",&
            description="Specifies a list of atoms belonging to the region.",&
            usage="LIST {integer} {integer} .. {integer}", required=.TRUE.,repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(region_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword=keyword,name="temperature",&
            description="The temperature in K used to "//&
            " initialize the velocities of the atoms in this region ",&
            usage="temperature 5.0",&
            default_r_val=cp_unit_to_cp2k(value=0.0_dp,unit_str="K",error=error),&
            unit_str="K", error=error)
       CALL section_add_keyword(region_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="temp_tol",&
            description="Maximum accepted  temperature deviation "//&
            " from the expected value for this region."//&
            "If  temp_tol=0 no rescaling is performed", &
            usage="temp_tol 0.0", default_r_val=0.0_dp, unit_str='K', error=error)
       CALL section_add_keyword(region_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section, region_section, error=error)
       CALL section_release(region_section,error=error)

       NULLIFY(print_key)
       CALL section_create(subsection,name="PRINT",&
            description="Collects all print_keys for thermal_regions",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL cp_print_key_section_create(print_key,"TEMPERATURE",&
            description="Controls output of temperature per region.", &
            print_level=high_print_level, common_iter_levels=1,&
            filename="",error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_thermal_region_section

! *****************************************************************************
!> \brief   Defines the parameters for the setup of a cascade simulation
!> \param   error variable to control error logging, stopping,...
!>          see module cp_error_handling
!> \author  Matthias Krack (MK)
!> \date    03.02.2012
!> \version 1.0
! *****************************************************************************
  SUBROUTINE create_cascade_section(section,error)

    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_cascade_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure = .FALSE.
    NULLIFY (keyword)
    NULLIFY (subsection)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN

       CALL section_create(section,name="CASCADE",&
                           description="Defines the parameters for the setup of a cascade simulation.",&
                           n_keywords=1,&
                           n_subsections=1,&
                           repeats=.FALSE.,&
                           required=.FALSE.,&
                           error=error)

       CALL keyword_create(keyword,&
                           name="_SECTION_PARAMETERS_",&
                           description="Controls the activation of the CASCADE section.",&
                           usage="&CASCADE on",&
                           default_l_val=.FALSE.,&
                           lone_keyword_l_val=.TRUE.,&
                           error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ENERGY",&
                           description="Total energy transferred to the system during the cascade event.",&
                           usage="ENERGY 20.0",&
                           default_r_val=0.0_dp,&
                           unit_str="keV",&
                           error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_create(subsection,name="ATOM_LIST",&
                           description="Defines a list of atoms for which the initial velocities are modified",&
                           n_keywords=1,&
                           n_subsections=0,&
                           repeats=.FALSE.,&
                           required=.TRUE.,&
                           error=error)

       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
                           description="Defines the list of atoms for which the velocities are modified. "//&
                                       "Each record consists of the atomic index, the velocity vector, and "//&
                                       "a weight to define which fraction of the total energy is assigned "//&
                                       "to the current atom.",&
                           usage="Atomic index   v(x) v(y) v(z)   weight",&
                           repeats=.TRUE.,&
                           required=.FALSE.,&
                           type_of_var=lchar_t,&
                           error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_cascade_section



END MODULE input_cp2k_motion
