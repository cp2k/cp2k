!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      - taken out of input_cp2k_motion
!> \author Ole Schuett
! *****************************************************************************

MODULE input_cp2k_md
  USE bibliography,                    ONLY: &
       Evans1983, Guidon2008, Kantorovich2008, Kantorovich2008a, Kolafa2004, &
       Kuhne2007, Minary2003, Ricci2003, Tuckerman1992, VandeVondele2002
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE f77_blas
  USE input_constants,                 ONLY: &
       add_last_numeric, debug_print_level, high_print_level, &
       isokin_ensemble, langevin_ensemble, low_print_level, &
       medium_print_level, npe_f_ensemble, npe_i_ensemble, nph_ensemble, &
       nph_uniaxial_damped_ensemble, nph_uniaxial_ensemble, npt_f_ensemble, &
       npt_i_ensemble, nve_ensemble, nvt_adiabatic_ensemble, nvt_ensemble, &
       reftraj_ensemble
  USE input_cp2k_barostats,            ONLY: create_barostat_section
  USE input_cp2k_thermostats,          ONLY: create_region_section,&
                                             create_thermo_fast_section,&
                                             create_thermo_slow_section,&
                                             create_thermostat_section
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
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_md'

PUBLIC :: create_md_section

CONTAINS

! *****************************************************************************
!> \brief ...
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
!> \brief Defines LANGEVIN section
!> \param section ...
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
            description="Controls the set of parameters to run a Langevin MD. "//&
            "The integrator used follows that given in the article by Ricci et al. "//&
            "The user can define regions in the system where the atoms inside "//&
            "undergoes Langevin MD, while those outside the regions undergoes "//&
            "NVE Born Oppenheimer MD. To define the regions, the user should "//&
            "use THERMAL_REGION subsection of MOTION%MD. The theory for "//&
            "Langevin MD involving sub-regions can be found in articles by "//&
            "Kantorovitch et al. All the references can be found in the links below.",&
            citations=(/Ricci2003, Kantorovich2008, Kantorovich2008a/),&
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
!> \param section ...
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
!> \brief ...
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
!> \brief section will contain some parameters for the shells dynamics
!> \param section ...
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
!> \brief section will contain some parameters for the adiabatic dynamics
!> \param section ...
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
!> \param section ...
!> \param error ...
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
         //"N = v / |v|;  y = x + delta * N;  F_t = F(y) - &lang; F(y) | N &rang; * N; "&
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
!> \brief input section used to define regions with different temperature
!>        initialization and control
!> \param section ...
!> \param error ...
!> \par History
!>   - Added input for langevin regions in thermal regions section
!>     (2014/02/04, LT)
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

       CALL section_create(section,name="thermal_region", &
            description="Define regions where different initialization and control "//&
            "of the temperature is used. When MOTION%MD%ENSEMBLE is set to LANGEVIN, "//&
            "this section controls if the atoms defined inside and outside the "//&
            "thermal regions should undergo Langevin MD or NVE Born-Oppenheimer MD. "//&
            "The theory behind Langevin MD using different regions can be found in "//&
            "articles by Kantorovitch et al. listed below.", &
            citations=(/Kantorovich2008, Kantorovich2008a/), &
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE., error=error)

       NULLIFY(region_section)
       NULLIFY(keyword,subsection)

       CALL keyword_create(keyword=keyword,name="force_rescaling", &
            description="Control the rescaling ot the velocities in all the regions, "//&
            "according to the temperature assigned to each reagion, when "//&
            "RESTART_VELOCITY in EXT_RESTART is active.", &
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword, name="do_langevin_default", &
            description="If ENSEMBLE is set to LANGEVIN, controls whether the "//&
            "atoms NOT defined in the thermal regions to undergo langevin MD "//&
            "or not. If not, then the atoms will undergo NVE Born-Oppenheimer MD.", &
            usage="do_langevin_default .FALSE.", &
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL section_create(region_section,name="DEFINE_REGION",&
            description="This section provides the possibility to define arbitrary region ", &
            n_keywords=3, n_subsections=0, repeats=.TRUE., required=.TRUE., error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword, name="LIST",&
            description="Specifies a list of atoms belonging to the region.",&
            usage="LIST {integer} {integer} .. {integer}", &
            required=.TRUE., repeats=.TRUE., n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(region_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,name="temperature",&
            description="The temperature in K used to initialize the velocities "//&
            "of the atoms in this region ", &
            usage="temperature 5.0", &
            default_r_val=cp_unit_to_cp2k(value=0.0_dp,unit_str="K",error=error), &
            unit_str="K", error=error)
       CALL section_add_keyword(region_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="temp_tol",&
            description="Maximum accepted temperature deviation from the expected "//&
            "value for this region. If temp_tol=0 no rescaling is performed", &
            usage="temp_tol 0.0", &
            default_r_val=0.0_dp, unit_str='K', error=error)
       CALL section_add_keyword(region_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="do_langevin", &
            description="When ENSEMBLE is set to LANGEVIN, Controls whether "//&
            "the atoms in the thermal region should undergo Langevin MD. If "//&
            "not, then they will undergo NVE Born-Oppenheimer MD.", &
            usage="do_langevin .TRUE.",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(region_section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL section_add_subsection(section, region_section, error=error)
       CALL section_release(region_section,error=error)

       NULLIFY(print_key)
       CALL section_create(subsection,name="PRINT",&
            description="Collects all print_keys for thermal_regions",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            required=.TRUE., error=error)

       CALL cp_print_key_section_create(print_key,"TEMPERATURE",&
            description="Controls output of temperature per region.", &
            print_level=high_print_level, common_iter_levels=1, &
            filename="", error=error)
       CALL section_add_subsection(subsection,print_key, error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"LANGEVIN_REGIONS",&
            description="Controls output of information on which atoms "//&
            "underwent Langevin MD and which atoms did not.", &
            print_level=high_print_level, &
            filename="", error=error)
       CALL section_add_subsection(subsection, print_key, error=error)
       CALL section_release(print_key, error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_thermal_region_section


! *****************************************************************************
!> \brief   Defines the parameters for the setup of a cascade simulation
!> \param section ...
!> \param   error variable to control error logging, stopping,...
!>          see module cp_error_handling
!> \date    03.02.2012
!> \author  Matthias Krack (MK)
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
                                       "to the current atom:"//&
                                       "<p><tt><big>Atomic_index v<sub>x</sub> v<sub>y</sub> v<sub>x</sub> "//&
                                       "Weight</big></tt></p>",&
                           usage="{{Integer} {Real} {Real} {Real} {Real}}",&
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


! *****************************************************************************
!> \brief Defines AVERAGES section
!> \param section ...
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
!> \param section ...
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

END MODULE input_cp2k_md
