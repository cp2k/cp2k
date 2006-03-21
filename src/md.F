!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
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
  USE cp_output_handling,              ONLY: cp_print_key_generate_filename
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: &
       Harris_MC_method, Harris_MD_method, dump_xmol, ensemble_name_types, &
       hmc_run, isokin_ensemble, langevin_ensemble, ldmc_run, &
       method_name_types, nph_ensemble, nph_uniaxial_damped_ensemble, &
       nph_uniaxial_ensemble, npt_f_ensemble, npt_i_ensemble, nve_ensemble, &
       nvt_ensemble, output_name_types
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE message_passing,                 ONLY: mp_bcast
  USE parallel_rng_types,              ONLY: next_random_number
  USE particle_types,                  ONLY: particle_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_md_section, simulation_parameters_type, &
       initialize_velocities, qs_md_parameters_type, simpar_units_type


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
     REAL (KIND=dp) :: tau_nhc
     REAL (KIND=dp) :: tau_cell
     REAL (KIND=dp), POINTER, DIMENSION ( : ) :: dt_yosh
     REAL (KIND=dp) :: shake_tol
     TYPE ( qs_md_parameters_type ) :: qs_md_param
     TYPE ( simpar_units_type ) :: units
     !MK Langevin
     REAL (KIND=dp) :: gamma
     REAL (KIND=dp) :: var_w
     !TdK HMC
     INTEGER        :: ld_steps
     REAL(KIND=dp)  :: temp_alpha
     REAL(KIND=dp)  :: temp_beta
     LOGICAL        :: tot_energy
     LOGICAL        :: ext_temp
     LOGICAL        :: qs_accept
     LOGICAL        :: rnd_velocities
     LOGICAL        :: vel_scaling
     LOGICAL        :: qs_start
     LOGICAL        :: restore_history
     LOGICAL        :: restore_momenta
     LOGICAL        :: semi_hybrid
     LOGICAL        :: soft_start
     LOGICAL        :: fill_wf_history
     LOGICAL        :: mmc_extrapolation
     LOGICAL        :: always_accept
     LOGICAL        :: drld
     LOGICAL        :: multiple_scf
     INTEGER        :: hmc_max_scf
     LOGICAL        :: lpsc
     INTEGER        :: md_method
     INTEGER        :: mc_method
  END TYPE simulation_parameters_type

  ! *** Global parameters ***
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'md'

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

  SUBROUTINE read_md_section ( simpar, globenv, error)
    TYPE(simulation_parameters_type), &
      INTENT(OUT)                            :: simpar
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER                                  :: iw
    INTEGER, DIMENSION(:), POINTER           :: iprint
    REAL(kind=dp)                            :: tmp_r
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: md_section, motion_section, &
                                                print_key
    TYPE(unit_convert_type), POINTER         :: unit

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    NULLIFY(unit, md_section, motion_section, print_key)
    CALL set_md_defaults(simpar, globenv)
    iw = globenv % scr

    motion_section => section_vals_get_subs_vals(globenv%input_file,&
         "MOTION",error=error)
    md_section => section_vals_get_subs_vals(motion_section,&
         "MD",error=error)
    CALL read_md_new_section(simpar, globenv,  md_section, error)

    IF (globenv%ionode.AND.(globenv%print_level >= 0)) THEN
       WRITE (iw, *)
    END IF

    IF (globenv%run_type_id == ldmc_run) THEN
       simpar%ensemble = langevin_ensemble
       simpar%rnd_velocities = .FALSE.
       simpar%vel_scaling = .FALSE.
       simpar%restore_history = .TRUE.
       simpar%restore_momenta = .FALSE.
       simpar%semi_hybrid = .TRUE.
    END IF

    IF (globenv%run_type_id == hmc_run .OR. globenv%run_type_id == ldmc_run) THEN

       IF (globenv%ionode.AND.(globenv%print_level >= 0)) THEN
          WRITE ( iw, '( A,T71,I10 )' ) ' HMC| Number of Harris LD steps per MC cycle ', &
               simpar%ld_steps
          WRITE ( iw, '( A,T71,F10.2 )' ) ' HMC| Temperature correction value alpha ', &
               simpar%temp_alpha
          WRITE ( iw, '( A,T71,F10.2 )' ) ' HMC| Temperature correction value beta ', &
               simpar%temp_beta
          WRITE ( iw, '( A,T75,L6,A )' ) ' HMC| TOT_ENERGY ', &
               simpar%tot_energy
          WRITE ( iw, '( A,T75,L6,A )' ) ' HMC| EXT_TEMP ', &
               simpar%ext_temp
          WRITE ( iw, '( A,T75,L6,A )' ) ' HMC| QS_ACCEPT ', &
               simpar%qs_accept
          WRITE ( iw, '( A,T75,L6,A )' ) ' HMC| RND_VELOCITIES ', &
               simpar%rnd_velocities
          WRITE ( iw, '( A,T75,L6,A )' ) ' HMC| VEL_SCALING ', &
               simpar%vel_scaling
          WRITE ( iw, '( A,T75,L6,A )' ) ' HMC| QS_START ', &
               simpar%qs_start
          WRITE ( iw, '( A,T75,L6,A )' ) ' HMC| RESTORE_HISTORY ', &
               simpar%restore_history
          WRITE ( iw, '( A,T75,L6,A )' ) ' HMC| RECOVER_MOMENTA ', &
               simpar%restore_momenta
          WRITE ( iw, '( A,T75,L6,A )' ) ' HMC| SEMI_HYBRID ', &
               simpar%semi_hybrid
          WRITE ( iw, '( A,T75,L6,A )' ) ' HMC| SOFT_START ', &
               simpar%soft_start
          WRITE ( iw, '( A,T75,L6,A )' ) ' HMC| FILL_WF_HISTORY ', &
               simpar%fill_wf_history
          WRITE ( iw, '( A,T75,L6,A )' ) ' HMC| MMC_EXTRAPOLATION ', &
               simpar%mmc_extrapolation
          WRITE ( iw, '( A,T75,L6,A )' ) ' HMC| ALWAYS_ACCEPT ', &
               simpar%always_accept
          WRITE ( iw, '( A,T75,L6,A )' ) ' HMC| DRLD ', &
               simpar%drld
          WRITE ( iw, '( A,T75,L6,A )' ) ' HMC| MULTIPLE_SCF ', &
               simpar%multiple_scf
          WRITE ( iw, '( A,T71,I10 )' ) ' HMC| HMC_MAX_SCF ', &
               simpar%hmc_max_scf
          WRITE ( iw, '( A,T75,L6,A )' ) ' HMC| LPSC ', &
               simpar%lpsc

          WRITE ( iw, '( A,T61,A )' ) ' HMC| MD Method ', &
               ADJUSTR ( method_name_types(simpar % md_method)(1:20))
          WRITE ( iw, '( A,T61,A )' ) ' HMC| MC Method ', &
               ADJUSTR ( method_name_types(simpar % mc_method)(1:20))

          WRITE ( iw, '( A,T20,A,T71,A )' ) ' HMC| File type', &
               'Print frequency[steps]', 'File names'
          print_key => section_vals_get_subs_vals(motion_section,"MD%PRINT%HMC",error=error)
          CALL section_vals_val_get(print_key,"EACH",i_vals=iprint,error=error)
          filename = cp_print_key_generate_filename(logger,print_key,&
               extension=".hmc",middle_name="pos",my_local=.FALSE.,error=error)          
          WRITE ( iw, '( A,T25,i5,T31,A50 )' ) ' HMC| HMC Statistics', iprint, &
               ADJUSTR ( TRIM(filename) )
          WRITE (iw, *)
       END IF
    END IF

    !MK Begin setup Langevin dynamics

    IF (simpar%ensemble == langevin_ensemble) THEN

      ! Normalization factor using a normal Gaussian random number distribution

      simpar%var_w = 2.0_dp*simpar%temp_ext*simpar%dt*simpar%gamma

      IF (globenv%ionode.AND.(globenv%print_level >= 0)) THEN
        tmp_r = cp_unit_from_cp2k(simpar%gamma,"fs^-1",error=error)
        WRITE (UNIT=iw,FMT="(T2,A,T71,ES10.3)")&
          "LD| Gamma [1/fs] ",tmp_r,&
          "LD| Variance [a.u.]",simpar%var_w
      END IF

    END IF

    !MK End setup Langevin dynamics

    !..write some information to output
    IF (globenv % ionode .AND. globenv % print_level>=0) THEN
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
               ' MD| Pressure [','Bar','] ', simpar % p_ext/unit%bar
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
       WRITE ( iw, '( A,T25,i5,T31,A50 )' ) ' MD| Coordinates',iprint(SIZE(iprint)), &
            ADJUSTR ( TRIM(filename) )

       print_key => section_vals_get_subs_vals(motion_section,"PRINT%VELOCITIES",error=error)
       NULLIFY(iprint)
       CALL section_vals_val_get(print_key,"EACH",i_vals=iprint,error=error)
       filename = cp_print_key_generate_filename(logger,print_key,&
            extension=".xyz",middle_name="vel",my_local=.FALSE.,error=error)
       WRITE ( iw, '( A,T25,i5,T31,A50 )' ) ' MD| Velocities',iprint(SIZE(iprint)), &
            ADJUSTR ( TRIM(filename) )

       print_key => section_vals_get_subs_vals(motion_section,"MD%PRINT%ENERGY",error=error)
       NULLIFY(iprint)
       CALL section_vals_val_get(print_key,"EACH",i_vals=iprint,error=error)
       filename =cp_print_key_generate_filename(logger,print_key,&
            extension=".ener",my_local=.FALSE.,error=error)
       WRITE ( iw, '( A,T25,i5,T31,A50 )' ) ' MD| Energies',iprint(SIZE(iprint)), &
            ADJUSTR ( TRIM(filename) )

       print_key => section_vals_get_subs_vals(motion_section,"PRINT%RESTART",error=error)
       NULLIFY(iprint)
       CALL section_vals_val_get(print_key,"EACH",i_vals=iprint,error=error)
       filename = cp_print_key_generate_filename(logger,print_key,&
            extension=".restart",my_local=.FALSE.,error=error)
       WRITE ( iw, '( A,T25,i5,T31,A50 )' ) ' MD| Dump',iprint(SIZE(iprint)), &
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
  SUBROUTINE read_md_new_section(simpar, globenv, md_section, error)
    TYPE(simulation_parameters_type), &
      INTENT(INOUT)                          :: simpar
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: md_section
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
       END IF

       tmp_section => section_vals_get_subs_vals(globenv%input_file,"GLOBAL",error=error)
       CALL section_vals_val_get(globenv%input_file,"MOTION%PRINT%TRAJECTORY%FORMAT", &
            i_val=simpar%output_format,error=error)

       tmp_section => section_vals_get_subs_vals(globenv%input_file,&
            "FORCE_EVAL%SUBSYS%CONSTRAINT",error=error)
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

       tmp_section => section_vals_get_subs_vals(md_section, "HMC", error=error)
       CALL section_vals_get(tmp_section, explicit=explicit, error=error)
       IF (explicit) THEN
          !CALL section_vals_val_get(tmp_section, "GAMMA", &
          !                          r_val=simpar%gamma, error=error)
          CALL section_vals_val_get(tmp_section, "LD_STEPS", &
               i_val=simpar%ld_steps, error=error)
          CALL section_vals_val_get(tmp_section, "TEMP_ALPHA", &
               r_val=simpar%temp_alpha, error=error)
          CALL section_vals_val_get(tmp_section, "TEMP_BETA", &
               r_val=simpar%temp_beta, error=error)
          CALL section_vals_val_get(tmp_section, "TOT_ENERGY", &
               l_val=simpar%tot_energy, error=error)
          CALL section_vals_val_get(tmp_section, "EXT_TEMP", &
               l_val=simpar%ext_temp, error=error)
          CALL section_vals_val_get(tmp_section, "QS_ACCEPT", &
               l_val=simpar%qs_accept, error=error)
          CALL section_vals_val_get(tmp_section, "RND_VELOCITIES", &
               l_val=simpar%rnd_velocities, error=error)
          CALL section_vals_val_get(tmp_section, "VEL_SCALING", &
               l_val=simpar%vel_scaling, error=error)
          CALL section_vals_val_get(tmp_section, "QS_START", &
               l_val=simpar%qs_start, error=error)
          CALL section_vals_val_get(tmp_section, "RESTORE_HISTORY", &
               l_val=simpar%restore_history, error=error)
          CALL section_vals_val_get(tmp_section, "RECOVER_MOMENTA", &
               l_val=simpar%restore_momenta, error=error)
          CALL section_vals_val_get(tmp_section, "SEMI_HYBRID", &
               l_val=simpar%semi_hybrid, error=error)
          CALL section_vals_val_get(tmp_section, "SOFT_START", &
               l_val=simpar%soft_start, error=error)
          CALL section_vals_val_get(tmp_section, "FILL_WF_HISTORY", &
               l_val=simpar%fill_wf_history, error=error)
          CALL section_vals_val_get(tmp_section, "MMC_extrapolation", &
               l_val=simpar%MMC_extrapolation, error=error)
          CALL section_vals_val_get(tmp_section, "ALWAYS_ACCEPT", &
               l_val=simpar%always_accept, error=error)
          CALL section_vals_val_get(tmp_section, "DRLD", &
               l_val=simpar%drld, error=error)
          CALL section_vals_val_get(tmp_section, "MULTIPLE_SCF", &
               l_val=simpar%multiple_scf, error=error)
          CALL section_vals_val_get(tmp_section, "HMC_MAX_SCF", &
               i_val=simpar%hmc_max_scf, error=error)
          CALL section_vals_val_get(tmp_section, "LPSC", &
               l_val=simpar%lpsc, error=error)

          CALL section_vals_val_get(tmp_section, "MD_METHOD", &
               i_val=simpar%md_method, error=error)
          CALL section_vals_val_get(tmp_section, "MC_METHOD", &
               i_val=simpar%mc_method, error=error)
       END IF

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
  SUBROUTINE set_md_defaults(simpar, globenv)
    TYPE(simulation_parameters_type), &
      INTENT(OUT)                            :: simpar
    TYPE(global_environment_type), POINTER   :: globenv

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
    simpar%var_w = 0.0_dp
    ! TdK HMC
    simpar%ld_steps = 10
    simpar%temp_alpha = 0.0_dp
    simpar%temp_beta = 1.0_dp
    simpar%tot_energy = .TRUE.
    simpar%ext_temp = .TRUE.
    simpar%qs_accept = .FALSE.
    simpar%rnd_velocities = .TRUE.
    simpar%vel_scaling = .FALSE.
    simpar%qs_start = .TRUE.
    simpar%restore_history = .TRUE.
    simpar%restore_momenta = .FALSE.
    simpar%semi_hybrid = .FALSE.
    simpar%soft_start = .FALSE.
    simpar%fill_wf_history = .TRUE.
    simpar%MMC_extrapolation = .FALSE.
    simpar%always_accept = .FALSE.
    simpar%drld = .FALSE.
    simpar%multiple_scf = .FALSE.
    simpar%hmc_max_scf = 1000
    simpar%lpsc = .FALSE.
    simpar%md_method = Harris_MD_method
    simpar%mc_method = Harris_MC_method

  END SUBROUTINE set_md_defaults





  !******************************************************************************
  !!****** md/initialize_velocities [1.0] *
  !!
  !!   NAME
  !!     initialize_velocities
  !!
  !!   SYNOPSIS
  !!     Subroutine initialize_velocities(simpar, part, globenv)
  !!       Implicit None
  !!       Type(simulation_parameters_type), Intent (IN):: simpar
  !!       Type(particle_type), Intent (INOUT), Dimension(:):: part
  !!       Type(global_environment_type), pointer:: globenv
  !!     End Subroutine initialize_velocities
  !!
  !!   FUNCTION
  !!     Initializes the velocities to the Maxwellian distribution
  !!
  !!   AUTHOR
  !!     CJM
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!*** *************************************************************************

  SUBROUTINE initialize_velocities ( simpar, part, globenv, work_section, error )


    TYPE(simulation_parameters_type), &
      INTENT(INOUT)                          :: simpar
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: part
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: work_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: i, iw, natoms
    LOGICAL                                  :: explicit
    REAL(KIND=dp)                            :: akin, denom, mass, mass_tot, &
                                                sc, temp
    REAL(KIND=dp), DIMENSION(3)              :: v, vcom, vcom_new
    REAL(KIND=dp), DIMENSION(:), POINTER     :: vel
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(unit_convert_type), POINTER         :: unit

!------------------------------------------------------------------------------

    iw = globenv % scr
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

    CALL get_cp2k_units ( unit ) 
    IF (ASSOCIATED(work_section)) &
         CALL section_vals_get(work_section,explicit=explicit,error=error)
    
    IF (explicit) THEN
       akin = 0.0_dp
       DO i = 1, natoms
          atomic_kind => part ( i ) % atomic_kind
          CALL get_atomic_kind ( atomic_kind = atomic_kind, mass = mass )
          mass_tot = mass_tot + mass
          IF (part(i)%is_fixed) THEN
             part ( i ) % v ( 1 ) = 0.0_dp
             part ( i ) % v ( 2 ) = 0.0_dp
             part ( i ) % v ( 3 ) = 0.0_dp
          ELSE
             CALL section_vals_val_get(work_section,"_DEFAULT_KEYWORD_",i_rep_val=i,r_vals=vel,error=error)
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
          IF (.NOT.part(i)%is_fixed) THEN
             v = part ( i ) % v ( : ) - vcom ( : )
             ! now compute the temp
             akin = akin + 0.5_dp * mass * SUM( v * v)
          END IF
       END DO
       temp = 0.0_dp
       IF (simpar % nfree /= 0) THEN 
          temp = 2.0_dp * akin / REAL ( simpar % nfree,KIND=dp)
       END IF
       vcom_new=vcom
    ELSE
       ! initializing velocities
       IF (globenv % ionode) THEN
          ! computing the COM velocities
          DO i = 1, natoms
             atomic_kind => part ( i ) % atomic_kind
             CALL get_atomic_kind ( atomic_kind = atomic_kind, mass = mass )
             mass_tot = mass_tot + mass
             IF (part(i)%is_fixed) THEN
                part ( i ) % v ( 1 ) = 0.0_dp
                part ( i ) % v ( 2 ) = 0.0_dp
                part ( i ) % v ( 3 ) = 0.0_dp
             ELSE
                part(i)%v(1) = next_random_number(globenv%gaussian_rng_stream)
                part(i)%v(2) = next_random_number(globenv%gaussian_rng_stream)
                part(i)%v(3) = next_random_number(globenv%gaussian_rng_stream)
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
             IF (.NOT.part(i)%is_fixed) THEN
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
             atomic_kind => part ( i ) % atomic_kind
             CALL get_atomic_kind ( atomic_kind = atomic_kind, mass = mass )
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
       CALL mp_bcast ( mass_tot, globenv % source, globenv % group )
       DO i = 1, natoms
          CALL mp_bcast ( part ( i ) % v, globenv % source, globenv % group )
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
  END SUBROUTINE initialize_velocities

END MODULE md
