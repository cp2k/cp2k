!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2006  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****** cp2k/md [1.0] *
!!
!!   NAME
!!     md
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE md
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp2k_units,                      ONLY: get_cp2k_units,&
                                             unit_convert_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_generate_filename,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE cp_linked_list_val,              ONLY: cp_sll_val_type, cp_sll_val_next
  USE input_val_types,                 ONLY: val_type, val_get
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: &
       dump_xmol, ensemble_name_types, isokin_ensemble, langevin_ensemble, &
       nph_ensemble, nph_uniaxial_damped_ensemble, nph_uniaxial_ensemble, &
       npt_f_ensemble, npt_i_ensemble, nve_ensemble, nvt_ensemble, output_name_types
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_list_get
  USE kinds,                           ONLY: default_path_length,&
                                             dp,&
                                             int_size
  USE message_passing,                 ONLY: mp_bcast
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE molecule_kind_types,             ONLY: fixd_constraint_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE parallel_rng_types,              ONLY: next_random_number
  USE particle_types,                  ONLY: particle_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'md'

  ! Public data types

  PUBLIC :: qs_md_parameters_type,&
            simpar_units_type,&
            simulation_parameters_type

  ! Public subroutines

  PUBLIC :: initialize_velocities,&
            read_md_section

  !!****s* md/qs_md_parameters_type [1.0] *
  !!
  !!   NAME
  !!     qs_md_parameters_type
  !!
  !!   FUNCTION
  !!     parameters for the quickstep md
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     - neighbor_lists_from_scratch: if the neighbor_lists should be build
  !!       from scratch
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!   MODIFICATION HISTORY
  !!     03.2003 created [fawzi]
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE qs_md_parameters_type
     LOGICAL :: neighbor_lists_from_scratch
  END TYPE qs_md_parameters_type
  !!***
  !****************************************************************************

  TYPE simpar_units_type
     CHARACTER ( LEN = 8 ) :: time, tau,            &
          verlet_skin, tau_cell,&
          press
  END TYPE simpar_units_type

  TYPE simulation_parameters_type
     INTEGER :: nsteps
     REAL (KIND=dp) :: dt
     REAL (KIND=dp) :: temp_ext
     REAL (KIND=dp) :: temp_baro
     REAL (KIND=dp) :: temp_tol
     REAL (KIND=dp) :: p_ext
     REAL (KIND=dp) :: cmass
     REAL (KIND=dp) :: v0
     REAL (KIND=dp) :: e0
     REAL (KIND=dp) :: v_shock
     REAL (KIND=dp) :: p0
     REAL (KIND=dp) :: f_annealing
     REAL (KIND=dp) :: gamma_nph
     INTEGER        :: ensemble
     INTEGER        :: output_format
     LOGICAL :: constraint
     LOGICAL :: annealing
     LOGICAL :: pv_numerical
     LOGICAL :: pv_availability
     INTEGER :: nfree
     INTEGER :: nc
     INTEGER :: nyosh
     INTEGER :: nhclen
     INTEGER :: info_constraint
     REAL (KIND=dp) :: tau_nhc
     REAL (KIND=dp) :: tau_cell
     REAL (KIND=dp), POINTER, DIMENSION ( : ) :: dt_yosh
     REAL (KIND=dp) :: shake_tol
     TYPE ( qs_md_parameters_type ) :: qs_md_param
     TYPE ( simpar_units_type )     :: units
     !MK Langevin
     REAL (KIND=dp) :: gamma
     REAL (KIND=dp) :: noisy_gamma
     REAL (KIND=dp) :: shadow_gamma
     REAL (KIND=dp) :: var_w
     LOGICAL        :: multi_time_switch, do_respa
     INTEGER        :: n_time_steps

  END TYPE simulation_parameters_type

  !!***
  !******************************************************************************

CONTAINS

  !******************************************************************************
  !!****** md/read_md_section [1.0] *
  !!
  !!   NAME
  !!     read_md_section
  !!
  !!   FUNCTION
  !!     reads the md section and prints out what has been read
  !!
  !!   AUTHOR
  !!     Teo, CJM
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!   NOTES
  !!*** ************************************************************************

  SUBROUTINE read_md_section ( simpar,root_section, para_env, globenv, subsys_section, error)
    TYPE(simulation_parameters_type), &
      INTENT(OUT)                            :: simpar
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER                                  :: iw
    INTEGER, DIMENSION(:), POINTER           :: iprint
    REAL(kind=dp)                            :: tmp_r1, tmp_r2, tmp_r3
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: md_section, motion_section, &
                                                print_key
    TYPE(unit_convert_type), POINTER         :: unit

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    NULLIFY(unit, md_section, motion_section, print_key)
    CALL set_md_defaults(simpar)

    motion_section => section_vals_get_subs_vals(root_section,&
         "MOTION",error=error)
    md_section => section_vals_get_subs_vals(motion_section,&
         "MD",error=error)

    logger => cp_error_get_logger(error)
    iw=cp_print_key_unit_nr(logger,md_section,"PRINT%PROGRAM_RUN_INFO",&
                                         extension=".log",error=error)

    CALL read_md_new_section(simpar, root_section, md_section, subsys_section, error)

    IF ( iw > 0 ) WRITE (iw, *)

    !MK Begin setup Langevin dynamics

    IF (simpar%ensemble == langevin_ensemble) THEN

      ! Normalization factor using a normal Gaussian random number distribution

      simpar%var_w = 2.0_dp*simpar%temp_ext*simpar%dt*(simpar%gamma+simpar%noisy_gamma)

      IF ( iw > 0 ) THEN
        tmp_r1 = cp_unit_from_cp2k(simpar%gamma,"fs^-1",error=error)
        tmp_r2 = cp_unit_from_cp2k(simpar%noisy_gamma, "fs^-1", error=error)
        tmp_r3 = cp_unit_from_cp2k(simpar%shadow_gamma, "fs^-1", error=error)
        WRITE (UNIT=iw,FMT="(T2,A,T71,ES10.3)")&
          "LD| Gamma [1/fs] ",tmp_r1,&
          "LD| Noisy Gamma [1/fs]",tmp_r2,&
          "LD| Shadow Gamma [1/fs]", tmp_r3,&
          "LD| Variance [a.u.]",simpar%var_w
      END IF

    END IF

    !MK End setup Langevin dynamics

    !..write some information to output
    IF ( iw > 0 ) THEN
       CALL get_cp2k_units ( unit )
       WRITE ( iw, '( A )' ) ' MD| Molecular Dynamics Protocol '
       WRITE ( iw, '( A,T61,A )' ) ' MD| Output Format ', &
            ADJUSTR ( output_name_types(simpar % output_format)(1:20))
       WRITE ( iw, '( A,T61,A )' ) ' MD| Ensemble Type ', &
            ADJUSTR ( ensemble_name_types(simpar % ensemble)(1:20))
       WRITE ( iw, '( A,T71,I10 )' ) ' MD| Number of Time Steps ', &
            simpar % nsteps
       WRITE ( iw, '( A,A2,A,T71,F10.2 )' ) ' MD| Time Step [','fs','] ',&
            simpar % dt /unit%fs
       WRITE ( iw, '( A,T71,F10.2 )' ) ' MD| Temperature [K] ', &
            simpar % temp_ext/ unit % kelvin
       WRITE ( iw, '( A,T71,F10.2 )' ) ' MD| Temperature tolerance [K] ', &
            simpar % temp_tol/ unit % kelvin
       IF ( simpar % annealing ) &
            WRITE ( iw, '( A,T71,F10.2 )' ) ' MD| Annealing ion factor      ', &
            simpar % f_annealing
       IF ( simpar % ensemble == nvt_ensemble .OR. simpar % ensemble == npt_i_ensemble &
            .OR. simpar % ensemble == npt_f_ensemble ) THEN
          WRITE ( iw, '( A )' ) ' MD| Nose-Hoover-Chain parameters '
          WRITE ( iw, '( A,T77,I4 )' ) ' MD| Nose-Hoover-Chain length ', &
               simpar % nhclen
          WRITE ( iw, '( A,A4,A,T71,F10.2 )' ) &
               ' MD| Nose-Hoover-Chain time constant [', 'fs' ,'] ',&
               simpar % tau_nhc/unit%fs
          WRITE ( iw, '( A,T77,I4 )' ) ' MD| Order of Yoshida integrator ', &
               simpar % nyosh
          WRITE ( iw, '( A,T77,I4 )' ) ' MD| Number of multiple time steps', &
               simpar % nc
       END IF
       IF (  simpar % ensemble == npt_i_ensemble .OR. simpar % ensemble == npt_f_ensemble ) THEN
          WRITE ( iw, '( A,A3,A, T71, F10.2 )' ) &
               ' MD| Pressure [','Bar','] ', simpar % p_ext*unit%bar
          WRITE ( iw, '( A,A4,A, T71, F10.2 )' ) &
               ' MD| Barostat time constant [','fs','] ', &
               simpar % tau_cell/unit%fs
       END IF
       IF (  simpar % ensemble == isokin_ensemble ) THEN
          WRITE ( iw, '( A )' ) ' MD| Simulation in the isokinetic ensemble'
       END IF
       IF ( simpar % constraint ) THEN
          WRITE ( iw, '( A )' ) ' MD| Constraints activated '
          WRITE ( iw, '( A,T71,G10.4 )' ) ' MD| Tolerance for shake ', &
               simpar % shake_tol
       END IF

       print_key => section_vals_get_subs_vals(motion_section,"MD%PRINT%PROGRAM_RUN_INFO",error=error)
       NULLIFY(iprint)
       CALL section_vals_val_get(print_key,"EACH",i_vals=iprint,error=error)
       WRITE ( iw, '( A,T63,i10,A )' ) ' MD| Print MD information every',iprint, ' step(s)'
       WRITE ( iw, '( A,T20,A,T71,A10 )' ) ' MD| File type','Print frequency[steps]', 'File names'

       print_key => section_vals_get_subs_vals(motion_section,"PRINT%TRAJECTORY",error=error)
       NULLIFY(iprint)
       CALL section_vals_val_get(print_key,"EACH",i_vals=iprint,error=error)
       filename = cp_print_key_generate_filename(logger,print_key,&
            extension=".xyz",middle_name="pos",my_local=.FALSE.,error=error)
       WRITE ( iw, '( A,T20,i10,T31,A50 )' ) ' MD| Coordinates',iprint(SIZE(iprint)), &
            ADJUSTR ( TRIM(filename) )

       IF ( simpar % ensemble == nph_ensemble .OR.&
               simpar % ensemble == nph_uniaxial_ensemble .OR.&
               simpar % ensemble == nph_uniaxial_damped_ensemble.OR.&
               simpar % ensemble == npt_i_ensemble.OR.&
               simpar % ensemble == npt_f_ensemble ) THEN

          print_key => section_vals_get_subs_vals(motion_section,"PRINT%CELL",error=error)
          NULLIFY(iprint)
          CALL section_vals_val_get(print_key,"EACH",i_vals=iprint,error=error)
          filename = cp_print_key_generate_filename(logger,print_key,&
              extension=".cell",my_local=.FALSE.,error=error)
         WRITE ( iw, '( A,T20,i10,T31,A50 )' ) ' MD| Simulation Cell',iprint(SIZE(iprint)), &
              ADJUSTR ( TRIM(filename) )
       END IF

       print_key => section_vals_get_subs_vals(motion_section,"PRINT%VELOCITIES",error=error)
       NULLIFY(iprint)
       CALL section_vals_val_get(print_key,"EACH",i_vals=iprint,error=error)
       filename = cp_print_key_generate_filename(logger,print_key,&
            extension=".xyz",middle_name="vel",my_local=.FALSE.,error=error)
       WRITE ( iw, '( A,T20,i10,T31,A50 )' ) ' MD| Velocities',iprint(SIZE(iprint)), &
            ADJUSTR ( TRIM(filename) )

       print_key => section_vals_get_subs_vals(motion_section,"MD%PRINT%ENERGY",error=error)
       NULLIFY(iprint)
       CALL section_vals_val_get(print_key,"EACH",i_vals=iprint,error=error)
       filename =cp_print_key_generate_filename(logger,print_key,&
            extension=".ener",my_local=.FALSE.,error=error)
       WRITE ( iw, '( A,T20,i10,T31,A50 )' ) ' MD| Energies',iprint(SIZE(iprint)), &
            ADJUSTR ( TRIM(filename) )

       print_key => section_vals_get_subs_vals(motion_section,"PRINT%RESTART",error=error)
       NULLIFY(iprint)
       CALL section_vals_val_get(print_key,"EACH",i_vals=iprint,error=error)
       filename = cp_print_key_generate_filename(logger,print_key,&
            extension=".restart",my_local=.FALSE.,error=error)
       WRITE ( iw, '( A,T20,i10,T31,A50 )' ) ' MD| Dump',iprint(SIZE(iprint)), &
            ADJUSTR ( TRIM(filename) )

       WRITE (iw,*)
       IF ( simpar % ensemble == nph_uniaxial_ensemble .OR.  &
            simpar % ensemble == nph_uniaxial_damped_ensemble ) THEN
          WRITE ( iw, '( A )' ) ' SHOCK| Uniaxial Shock Parameters: '
          WRITE ( iw, '( A,A4,A,T71,F10.4 )' ) &
               ' SHOCK| Shock Velocity [', 'm/s' ,'] ',&
               simpar % v_shock * 1.e-10_dp * 1.e15_dp / unit % angstrom * unit%fs
          WRITE ( iw, '( A,A4,A,T71,F10.4 )' ) &
               ' SHOCK| Damping Coefficient [', '1/fs' ,'] ',&
               simpar % gamma_nph*unit%fs
          WRITE ( iw, '( A,A3,A, T71, F10.2 )' ) &
               ' SHOCK| Pressure [','Bar','] ', simpar % p0*unit%bar
          WRITE ( iw, '( A,A4,A, T71, E10.4 )' ) &
               ' SHOCK| Barostat Mass [', 'au', '] ', &
               simpar % cmass
       END IF
    END IF

    CALL cp_print_key_finished_output(iw,logger,md_section,&
               "PRINT%PROGRAM_RUN_INFO", error=error)

  END SUBROUTINE read_md_section

  !******************************************************************************
  !!****** md/read_md_new_section [1.0] *
  !!
  !!   NAME
  !!     read_md_new_section
  !!
  !!
  !!   FUNCTION
  !!     Parses the MD input section
  !!
  !!   AUTHOR
  !!     teo
  !!
  !!   MODIFICATION HISTORY
  !!
  !!
  !!*** *************************************************************************
  SUBROUTINE read_md_new_section(simpar, root_section, md_section, subsys_section, error)
    TYPE(simulation_parameters_type), &
      INTENT(INOUT)                          :: simpar
    TYPE(section_vals_type), POINTER         :: root_section, md_section, &
                                                subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_md_new_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: explicit, failure
    TYPE(section_vals_type), POINTER         :: tmp_section

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CALL section_vals_val_get(md_section,"ANNEALING",r_val=simpar%f_annealing,error=error)
       IF (simpar%f_annealing /= 1.0_dp) simpar % annealing = .TRUE.
       CALL section_vals_val_get(md_section,"ENSEMBLE",i_val=simpar%ensemble,error=error)
       CALL section_vals_val_get(md_section,"STEPS",i_val=simpar%nsteps,error=error)
       CALL section_vals_val_get(md_section,"TIMESTEP",r_val=simpar%dt,error=error)
       CALL section_vals_val_get(md_section,"TEMPERATURE",r_val=simpar%temp_ext,error=error)
       CALL section_vals_val_get(md_section,"TEMP_TOL",r_val=simpar%temp_tol,error=error)
       CALL section_vals_val_get(md_section,"PV_AVA",l_val=simpar%pv_availability,error=error)
       CALL section_vals_val_get(md_section,"PV_NUM",l_val=simpar%pv_numerical,error=error)

       IF (simpar%ensemble == langevin_ensemble) THEN
          CALL section_vals_val_get(md_section, "GAMMA", r_val=simpar%gamma, error=error)
          CALL section_vals_val_get(md_section, "NOISY_GAMMA", r_val=simpar%noisy_gamma, &
                                    error=error)
          CALL section_vals_val_get(md_section, "SHADOW_GAMMA", r_val=simpar%shadow_gamma, &
                                    error=error)
       END IF

       tmp_section => section_vals_get_subs_vals(root_section,"GLOBAL",error=error)
       CALL section_vals_val_get(root_section,"MOTION%PRINT%TRAJECTORY%FORMAT", &
            i_val=simpar%output_format,error=error)

       tmp_section => section_vals_get_subs_vals(root_section,&
            "MOTION%CONSTRAINT",error=error)
       CALL section_vals_get(tmp_section,explicit=explicit,error=error)
       IF (explicit) THEN
          simpar%constraint = .TRUE.
          CALL section_vals_val_get(tmp_section,"SHAKE",r_val=simpar%shake_tol,error=error)
       END IF

       IF (simpar%ensemble == nph_uniaxial_ensemble.OR.simpar%ensemble == nph_uniaxial_damped_ensemble) THEN
          tmp_section => section_vals_get_subs_vals(md_section,"UNIAXIAL",error=error)
          CALL section_vals_val_get(tmp_section,"PRESSURE",r_val=simpar%p0,error=error)
          CALL section_vals_val_get(tmp_section,"ENERGY",r_val=simpar%e0,error=error)
          CALL section_vals_val_get(tmp_section,"VOLUME",r_val=simpar%v0,error=error)
          CALL section_vals_val_get(tmp_section,"GAMMA",r_val=simpar%gamma_nph,error=error)
          IF ( simpar % gamma_nph /= 0.0_dp ) simpar % ensemble = nph_uniaxial_damped_ensemble
          CALL section_vals_val_get(tmp_section,"CMASS",r_val=simpar%cmass,error=error)
          CALL section_vals_val_get(tmp_section,"VSHOCK",r_val=simpar%v_shock,error=error)
       END IF

       tmp_section => section_vals_get_subs_vals(md_section,"BAROSTAT",error=error)
       CALL section_vals_val_get(tmp_section,"PRESSURE",r_val=simpar%p_ext,error=error)
       CALL section_vals_val_get(tmp_section,"TIMECON",r_val=simpar%tau_cell,error=error)

       tmp_section => section_vals_get_subs_vals(md_section,"NOSE",error=error)
       CALL section_vals_val_get(tmp_section,"LENGTH",i_val=simpar%nhclen,error=error)
       CALL section_vals_val_get(tmp_section,"YOSHIDA",i_val=simpar%nyosh,error=error)
       CALL section_vals_val_get(tmp_section,"TIMECON",r_val=simpar%tau_nhc,error=error)
       CALL section_vals_val_get(tmp_section,"MTS",i_val=simpar%nc,error=error)

       tmp_section => section_vals_get_subs_vals(md_section,"RESPA",error=error)
       CALL section_vals_get(tmp_section, explicit=simpar%do_respa, error=error)
       CALL section_vals_val_get(tmp_section,"FREQUENCY",i_val=simpar%n_time_steps,&
            ignore_required=.TRUE.,error=error)
       simpar%multi_time_switch = .FALSE.
       IF (simpar%do_respa) simpar%multi_time_switch = .TRUE.

       tmp_section => section_vals_get_subs_vals(md_section,"QS_MD",error=error)
       CALL section_vals_get(tmp_section,explicit=explicit,error=error)
       IF (explicit) THEN
          CALL section_vals_val_get(tmp_section,"NEIGHBOR_LISTS_FROM_SCRATCH",&
               l_val=simpar%qs_md_param%neighbor_lists_from_scratch, error=error)
       END IF
    END IF

  END SUBROUTINE read_md_new_section

  !******************************************************************************
  !!****** md/set_md_defaults [1.0] *
  !!
  !!   NAME
  !!     set_md_defaults
  !!
  !!
  !!   FUNCTION
  !!     Initializes the simpar
  !!
  !!   AUTHOR
  !!     teo
  !!
  !!   MODIFICATION HISTORY
  !!
  !!
  !!*** *************************************************************************
  SUBROUTINE set_md_defaults(simpar)
    TYPE(simulation_parameters_type), &
      INTENT(OUT)                            :: simpar

    simpar % units % time = 'FS'
    simpar % units % press = 'BAR'
    simpar % units % tau = 'FS'
    simpar % units % tau_cell = 'FS'
    simpar % units % verlet_skin = 'ANGSTROMS'
    simpar % output_format = dump_xmol
    simpar % ensemble = nve_ensemble
    simpar % nsteps = 0
    simpar % dt = 0.0_dp
    simpar % temp_ext = 0.0_dp
    simpar % temp_baro = 0.0_dp
    simpar % temp_tol = 0.0_dp
    simpar % p_ext = 0.0_dp
    simpar % f_annealing = 1.0_dp
    simpar % constraint = .FALSE.
    simpar % annealing = .FALSE.
    simpar % pv_numerical=.FALSE.
    simpar % pv_availability=.FALSE.
    simpar % shake_tol = 1.0E-6_dp
    simpar % nhclen = 3
    simpar % nc = 2
    simpar % nyosh = 3
    simpar % tau_nhc = 1000.0_dp
    simpar % tau_cell = 1000.0_dp
    simpar % qs_md_param % neighbor_lists_from_scratch = .TRUE.
    simpar % p0 = 0._dp
    simpar % v0 = 0._dp
    simpar % e0 = 0._dp
    simpar % v_shock = 0.0_dp
    simpar % gamma_nph = 0.0_dp
    simpar % cmass = 0.0_dp
    !MK Langevin
    simpar%gamma = 0.0_dp
    simpar%noisy_gamma = 0.0_dp
    simpar%shadow_gamma = 0.0_dp
    simpar%var_w = 0.0_dp

  END SUBROUTINE set_md_defaults

  !******************************************************************************
  !!****** md/initialize_velocities [1.0] *
  !!
  !!   NAME
  !!     initialize_velocities
  !!
  !!   SYNOPSIS
  !!     Subroutine initialize_velocities(simpar, part)
  !!       Implicit None
  !!       Type(simulation_parameters_type), Intent (IN):: simpar
  !!       Type(particle_type), Intent (INOUT), Dimension(:):: part
  !!     End Subroutine initialize_velocities
  !!
  !!   FUNCTION
  !!     Initializes the velocities to the Maxwellian distribution
  !!
  !!   AUTHOR
  !!     CJM,MK
  !!
  !!   MODIFICATION HISTORY
  !!     - is_fixed removed from particle_type
  !!
  !!*** *************************************************************************

  SUBROUTINE initialize_velocities(simpar,part,para_env,globenv,molecule_kinds,&
                                   print_section,work_section,force_rescaling,error)

    TYPE(simulation_parameters_type), &
      INTENT(INOUT)                          :: simpar
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: part
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(section_vals_type), POINTER         :: print_section, work_section
    LOGICAL, INTENT(IN), OPTIONAL            :: force_rescaling
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_velocities'

    INTEGER                                  :: handle, i, ifixd, &
                                                imolecule_kind, istat, iw, &
                                                natoms, nfixd
    LOGICAL                                  :: explicit, is_ok
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: is_fixed
    REAL(KIND=dp)                            :: akin, denom, mass, mass_tot, &
                                                sc, temp
    REAL(KIND=dp), DIMENSION(3)              :: v, vcom, vcom_new
    REAL(KIND=dp), DIMENSION(:), POINTER     :: vel
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(unit_convert_type), POINTER         :: unit
    TYPE(cp_sll_val_type), POINTER           :: list
    TYPE(val_type), POINTER                  :: val

! -------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    logger => cp_error_get_logger(error)
    iw=cp_print_key_unit_nr(logger,print_section,"PROGRAM_RUN_INFO",&
                                         extension=".log",error=error)

    IF (iw>0) THEN
       WRITE ( iw, '( A, A, A )' ) ' **********************', &
            ' begin of velocity initialization ', '***********************'
    END IF

    ! Initializing parameters
    vcom = 0.0_dp
    mass_tot = 0.0_dp
    vcom_new = 0.0_dp
    denom = 0.0_dp
    natoms = SIZE ( part )
    explicit = .FALSE.

    ! Build a list of all fixed atoms (if there any)

    ALLOCATE (is_fixed(natoms),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "is_fixed",natoms*int_size)
    is_fixed(:) = .FALSE.

    molecule_kind_set => molecule_kinds%els

    DO imolecule_kind=1,molecule_kinds%n_els
      molecule_kind => molecule_kind_set(imolecule_kind)
      CALL get_molecule_kind(molecule_kind=molecule_kind,&
                             nfixd=nfixd,&
                             fixd_list=fixd_list)
      DO ifixd=1,nfixd
         IF (.NOT.fixd_list(ifixd)%restraint%active) is_fixed(fixd_list(ifixd)%fixd) = .TRUE.
      END DO
    END DO

    ! Read or initialize the particle velocities

    CALL get_cp2k_units ( unit )
    IF (ASSOCIATED(work_section)) &
         CALL section_vals_get(work_section,explicit=explicit,error=error)

    IF (explicit) THEN
       akin = 0.0_dp
       CALL section_vals_list_get(work_section,"_DEFAULT_KEYWORD_",list=list,error=error)
       DO i = 1, natoms
          is_ok=cp_sll_val_next(list,val,error=error)
          CALL val_get(val,r_vals=vel,error=error)
          atomic_kind => part ( i ) % atomic_kind
          CALL get_atomic_kind ( atomic_kind = atomic_kind, mass = mass )
          mass_tot = mass_tot + mass
          IF (is_fixed(i)) THEN
             part ( i ) % v ( 1 ) = 0.0_dp
             part ( i ) % v ( 2 ) = 0.0_dp
             part ( i ) % v ( 3 ) = 0.0_dp
          ELSE
             ! CALL section_vals_val_get(work_section,"_DEFAULT_KEYWORD_",i_rep_val=i,r_vals=vel,error=error)
             part ( i ) % v  = vel
             vcom  = vcom + mass * part ( i ) % v
             denom = denom + mass
          END IF
       END DO
       vcom = vcom / denom
       DO i = 1, natoms
          atomic_kind => part ( i ) % atomic_kind
          CALL get_atomic_kind ( atomic_kind = atomic_kind, mass = mass )
          ! first subtract the COM
          IF (.NOT.is_fixed(i)) THEN
             v = part ( i ) % v ( : ) - vcom ( : )
             ! now compute the temp
             akin = akin + 0.5_dp * mass * SUM( v * v)
          END IF
       END DO
       temp = 0.0_dp
       IF (simpar % nfree /= 0) THEN
          temp = 2.0_dp * akin / REAL ( simpar % nfree,KIND=dp)
          IF (PRESENT(force_rescaling)) THEN
            IF (force_rescaling) THEN
               sc = SQRT ( simpar % temp_ext / temp )
               DO i = 1, natoms
                  part ( i ) % v ( 1 ) = sc * part ( i ) % v ( 1 )
                  part ( i ) % v ( 2 ) = sc * part ( i ) % v ( 2 )
                  part ( i ) % v ( 3 ) = sc * part ( i ) % v ( 3 )
               END DO
            END IF
          END IF
       END IF
       vcom_new=vcom
    ELSE
       ! initializing velocities
       IF (para_env%ionode) THEN
          ! computing the COM velocities
          DO i = 1, natoms
             atomic_kind => part ( i ) % atomic_kind
             CALL get_atomic_kind ( atomic_kind = atomic_kind, mass = mass )
             mass_tot = mass_tot + mass
             IF (is_fixed(i)) THEN
                part ( i ) % v ( 1 ) = 0.0_dp
                part ( i ) % v ( 2 ) = 0.0_dp
                part ( i ) % v ( 3 ) = 0.0_dp
             ELSE
                part(i)%v(1) = next_random_number(globenv%gaussian_rng_stream,error=error)
                part(i)%v(2) = next_random_number(globenv%gaussian_rng_stream,error=error)
                part(i)%v(3) = next_random_number(globenv%gaussian_rng_stream,error=error)
                vcom ( 1 ) = vcom ( 1 ) + mass * part ( i ) % v ( 1 )
                vcom ( 2 ) = vcom ( 2 ) + mass * part ( i ) % v ( 2 )
                vcom ( 3 ) = vcom ( 3 ) + mass * part ( i ) % v ( 3 )
                denom = denom + mass
             END IF
          END DO
          vcom = vcom / denom

          akin = 0.0_dp
          DO i = 1, natoms
             atomic_kind => part ( i ) % atomic_kind
             CALL get_atomic_kind ( atomic_kind = atomic_kind, mass = mass )
             ! first subtract the COM
             IF (.NOT.is_fixed(i)) THEN
                part ( i ) % v ( : ) = part ( i ) % v ( : ) - vcom ( : )
                ! now compute the temp
                akin = akin + 0.5_dp *  &
                     ( mass * part ( i ) % v ( 1 ) * part ( i ) % v ( 1 ) )
                akin = akin + 0.5_dp *  &
                     ( mass * part ( i ) % v ( 2 ) * part ( i ) % v ( 2 ) )
                akin = akin + 0.5_dp *  &
                     ( mass * part ( i ) % v ( 3 ) * part ( i ) % v ( 3 ) )
             END IF
          END DO
          temp = 0.0_dp
          sc   = 0.0_dp
          IF (simpar % nfree /= 0) THEN
             temp = 2.0_dp * akin / REAL ( simpar % nfree,KIND=dp)
             ! scale velocities to get the correct initial temperature and
             sc = SQRT ( simpar % temp_ext / temp )
          END IF
          DO i = 1, natoms
             part ( i ) % v ( 1 ) = sc * part ( i ) % v ( 1 )
             part ( i ) % v ( 2 ) = sc * part ( i ) % v ( 2 )
             part ( i ) % v ( 3 ) = sc * part ( i ) % v ( 3 )
          END DO

          ! Moving the velocities so that the COM velocity is zero
          akin = 0.0_dp
          vcom_new = 0.0_dp
          DO i = 1, natoms
             atomic_kind => part ( i ) % atomic_kind
             CALL get_atomic_kind ( atomic_kind = atomic_kind, mass = mass )
             vcom_new ( : ) = vcom_new ( : ) + mass * part ( i ) % v ( : )
             akin = akin + 0.5_dp *  &
                  ( mass * part ( i ) % v ( 1 ) * part ( i ) % v ( 1 ) )
             akin = akin + 0.5_dp *  &
                  ( mass * part ( i ) % v ( 2 ) * part ( i ) % v ( 2 ) )
             akin = akin + 0.5_dp *  &
                  ( mass * part ( i ) % v ( 3 ) * part ( i ) % v ( 3 ) )
          END DO
       END IF
#if defined(__parallel)
       CALL mp_bcast ( mass_tot, para_env%source, para_env%group )
       DO i = 1, natoms
          CALL mp_bcast ( part ( i ) % v, para_env%source, para_env%group )
       END DO
#endif
    END IF

    IF (iw>0) THEN
       IF (simpar % nfree == 0) THEN
          temp = 2.0_dp * akin / REAL ( 1,KIND=dp) ! Should always be zero...
       ELSE
          temp = 2.0_dp * akin / REAL ( simpar % nfree,KIND=dp)
       END IF
       WRITE ( iw, '( A, T61, F18.2, A2 )' ) &
            ' Initial Temperature ', temp / unit % kelvin, " K"
       WRITE ( iw, '( A, T61, F20.12 )' ) &
            ' Centre of mass velocity in direction x:', vcom_new ( 1 )
       WRITE ( iw, '( A, T61, F20.12 )' ) &
            ' Centre of mass velocity in direction y:', vcom_new ( 2 )
       WRITE ( iw, '( A, T61, F20.12 )' ) &
            ' Centre of mass velocity in direction z:', vcom_new ( 3 )
       WRITE ( iw, '( A, A, A, / )' ) ' ***********************', &
            ' end of velocity initialization ', '************************'
    END IF

    IF ( simpar % ensemble == nph_ensemble .OR.&
         simpar % ensemble == nph_uniaxial_ensemble .OR.&
         simpar % ensemble == nph_uniaxial_damped_ensemble ) &
         simpar % v_shock = simpar % v_shock * SQRT ( mass_tot )

    DEALLOCATE (is_fixed,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"is_fixed")

    CALL cp_print_key_finished_output(iw,logger,print_section,&
               "PROGRAM_RUN_INFO", error=error)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE initialize_velocities

  ! ***************************************************************************

END MODULE md
