!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

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
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_logger_type,&
                                             cp_warning_level
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE cp2k_units,                      ONLY: get_cp2k_units,&
                                             unit_convert_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: nve_ensemble,&
                                             nvt_ensemble,&
                                             npt_i_ensemble,&
                                             npt_f_ensemble,&
                                             nph_ensemble,&
                                             nph_uniaxial_ensemble,&
                                             nph_uniaxial_damped_ensemble,&
                                             isokin_ensemble, &
                                             reftraj_ensemble,&
                                             langevin_ensemble,&
                                             nhc_global,&
                                             nhc_massive,&
                                             ensemble_name_types,&
                                             ensemble_name_index,&
                                             md_init_start,&
                                             md_pos_restart,&
                                             md_pos_vel_restart,&
                                             md_all_restart,&
                                             restart_name_types,&
                                             restart_name_index,&
                                             dump_atomic,&
                                             dump_xmol,&
                                             dump_dcd,&
                                             output_name_types,&
                                             output_name_index,&
                                             method_name_types, &
                                             hmc_run, &
                                             ldmc_run, &
                                             FIST_MD_method, &
                                             KG_MD_method, &
                                             QS_MD_method, &
                                             Harris_MD_method, &
                                             KG_MC_method, &
                                             QS_MC_method, &
                                             Harris_MC_method
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_get,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_bcast
  USE particle_types,                  ONLY: particle_type
  USE string_utilities,                ONLY: uppercase,&
                                             xstring
  USE util,                            ONLY: gasdev
  USE cp_units, only: cp_unit_from_cp2k
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_md_section, simulation_parameters_type, &
            initialize_velocities, mdio_parameters_type, &
            qs_md_parameters_type, simpar_units_type


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
  END TYPE 

  TYPE simulation_parameters_type
     INTEGER :: nsteps
     REAL (KIND=dp) :: walltime ! target maximum wall time of a simulation (alternative to nsteps)
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
     INTEGER        :: read_type
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
     INTEGER        :: md_method
     INTEGER        :: mc_method
  END TYPE simulation_parameters_type

  TYPE mdio_parameters_type
     INTEGER :: icrd, icoef, ivel, iptens, iener, itemp, idump, iscreen, ipot
     CHARACTER ( LEN = 50 ) :: crd_file_name
     CHARACTER ( LEN = 50 ) :: vel_file_name
     CHARACTER ( LEN = 50 ) :: ptens_file_name
     CHARACTER ( LEN = 50 ) :: pot_file_name
     CHARACTER ( LEN = 50 ) :: ener_file_name
     CHARACTER ( LEN = 50 ) :: temp_file_name
     CHARACTER ( LEN = 50 ) :: dump_file_name
     CHARACTER ( LEN = 50 ) :: coef_file_name
     CHARACTER ( LEN = 50 ) :: hmc_file_name
  END TYPE mdio_parameters_type

! *** Global parameters ***
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='md'

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

  SUBROUTINE read_md_section ( simpar, globenv, mdio, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(simulation_parameters_type), &
         INTENT(OUT)                            :: simpar
    TYPE(global_environment_type), &
         INTENT(IN)                             :: globenv
    TYPE(mdio_parameters_type), INTENT(OUT)     :: mdio
    TYPE(cp_error_type), INTENT(inout)          :: error
    ! Local Variables
    CHARACTER(LEN=20)                        :: string
    CHARACTER(LEN=5)                         :: label
    INTEGER                                  :: ia, ie, ierror, ILEN, iw
    TYPE(section_vals_type), POINTER         :: md_section
    TYPE(unit_convert_type), POINTER         :: unit
    REAL(kind=dp) :: tmp_r

    NULLIFY(unit, md_section)
    CALL set_md_defaults(simpar, globenv, mdio)
    iw = globenv % scr

    md_section => section_vals_get_subs_vals(globenv%input_file,&
       "MOTION%MD",error=error)
    CALL read_md_new_section(simpar, globenv,  md_section, error)

!TdK Turn on Langevin dynamics in case HMC is activated in the harris section
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

       WRITE ( iw, '( A,T61,A )' ) ' HMC| MD Method ', &
            ADJUSTR ( method_name_types(simpar % md_method)(1:20))
       WRITE ( iw, '( A,T61,A )' ) ' HMC| MC Method ', &
            ADJUSTR ( method_name_types(simpar % mc_method)(1:20))

       WRITE ( iw, '( A,T20,A,T71,A )' ) ' HMC| File type', &
            'Print frequency[steps]', 'File names'
       WRITE ( iw, '( A,T25,i5,T31,A )' ) ' HMC| HMC Statistics', 1, &
            ADJUSTR ( mdio % hmc_file_name ( 1:50 ) )
       WRITE (iw, *)
     END IF
  END IF

!MK begin setup Langevin
  IF (simpar%ensemble == langevin_ensemble) THEN
  
    IF (globenv%ionode.AND.(globenv%print_level >= 0)) THEN
       tmp_r=cp_unit_from_cp2k(simpar%gamma,"fs^-1",error=error)
      WRITE (iw,"(T2,A,T71,ES10.3)")&
        "LD| Gamma [1/fs] ",tmp_r
    END IF
  
    !MK Normalization factor uniform random number distribution [-0.5,0.5]
  
    simpar%var_w = 24.0_dp*simpar%temp_ext*simpar%dt*simpar%gamma
  
    !MK Normalization factor Gaussian random number distribution
    !MK simpar%var_w = 2.0_dp*simpar%temp_ext*simpar%dt*simpar%gamma
  
    IF (globenv%ionode.AND.(globenv%print_level >= 0)) THEN
      WRITE ( iw, '( A,T65,ES16.8 )' ) ' LD| External temperature T [a.u.] ', &
            simpar % temp_ext
      WRITE ( iw, '( A,T65,ES16.8 )' ) ' LD| Timestep dt [a.u.] ', &
            simpar % dt
      WRITE ( iw, '( A,T65,ES16.8 )' ) ' LD| Gamma [a.u.] ', &
            simpar % gamma
      !WRITE ( iw, '( A,T65,ES16.8 )' ) ' LD| W ', &
      !      simpar % var_w
      WRITE ( iw, *)
    END IF
  
  END IF
!MK End setup Langevin

!..write some information to output
    IF (globenv % ionode .AND. globenv % print_level>=0) THEN
       CALL get_cp2k_units ( unit )
       WRITE ( iw, '( A )' ) ' MD| Molecular Dynamics Protocol '
       WRITE ( iw, '( A,T61,A )' ) ' MD| Output Format ', &
            ADJUSTR ( output_name_types(simpar % output_format)(1:20))
       WRITE ( iw, '( A,T61,A )' ) ' MD| Ensemble Type ', & 
            ADJUSTR ( ensemble_name_types(simpar % ensemble)(1:20))
       WRITE ( iw, '( A,T61,A )' ) ' MD| Restart Type ', &
            ADJUSTR ( restart_name_types(simpar % read_type)(1:20))
       WRITE ( iw, '( A,T71,I10 )' ) ' MD| Number of Time Steps ', &
            simpar % nsteps
       WRITE ( iw, '( A,T71,F10.2 )' ) ' MD| Target Walltime [sec]', &
            simpar % walltime
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
          WRITE ( iw, '( A,A4,A,T71,F10.4 )' ) &
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
          WRITE ( iw, '( A,A4,A, T71, F10.4 )' ) &
               ' MD| Barostat time constant [','fs','] ', &
               simpar % tau_cell*unit%fs
       END IF
       IF (  simpar % ensemble == isokin_ensemble ) THEN
          WRITE ( iw, '( A )' ) ' MD| Simulation in the isokinetic ensemble'
       END IF
       IF ( simpar % constraint ) THEN
          WRITE ( iw, '( A )' ) ' MD| Constraints activated '
          WRITE ( iw, '( A,T71,G10.4 )' ) ' MD| Tolerance for shake ', &
               simpar % shake_tol
       END IF
       
       WRITE ( iw, '( A,T63,i10,A )' ) ' MD| Print MD information every', &
            mdio % iscreen, ' step(s)'
       WRITE ( iw, '( A,T20,A,T71,A )' ) ' MD| File type', &
            'Print frequency[steps]', 'File names'
       WRITE ( iw, '( A,T25,i5,T31,A )' ) ' MD| Coordinates', mdio % icrd, &
            ADJUSTR ( mdio % crd_file_name ( 1:50 ) )
       WRITE ( iw, '( A,T25,i5,T31,A )' ) ' MD| Velocities', mdio % ivel, &
            ADJUSTR ( mdio % vel_file_name ( 1:50 ) )
       WRITE ( iw, '( A,T25,i5,T31,A )' ) ' MD| Energies', mdio % iener, &
            ADJUSTR ( mdio % ener_file_name ( 1:50 ) )
       WRITE ( iw, '( A,T25,i5,T31,A )' ) ' MD| Potential Energies',     &
            mdio % ipot, ADJUSTR ( mdio % pot_file_name ( 1:50 ) )
       WRITE ( iw, '( A,T25,i5,T31,A )' ) ' MD| Temperature', mdio % itemp, &
            ADJUSTR ( mdio % temp_file_name ( 1:50 ) )
       WRITE ( iw, '( A,T25,i5,T31,A )' ) ' MD| Pressure Tensor', mdio % iptens, &
            ADJUSTR ( mdio % ptens_file_name ( 1:50 ) )
       WRITE ( iw, '( A,T25,i5,T31,A )' ) ' MD| Dump', mdio % idump, &
            ADJUSTR ( mdio % dump_file_name ( 1:50 ) )
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
  IMPLICIT NONE
  ! Arguments
  TYPE(simulation_parameters_type), &
       INTENT(INOUT)                          :: simpar
  TYPE(global_environment_type), &
       INTENT(IN)                             :: globenv
  TYPE(cp_error_type), INTENT(inout)          :: error
  TYPE(section_vals_type), POINTER            :: md_section
  ! Local Variables
  CHARACTER(len=*), PARAMETER :: routineN = 'read_md_section', &
       routineP = moduleN//':'//routineN
  TYPE(section_vals_type), POINTER            :: tmp_section
  LOGICAL :: explicit, failure
  INTEGER :: n_rep

  failure = .FALSE.
  IF (.NOT.failure) THEN
     CALL section_vals_val_get(md_section,"ANNEALING",r_val=simpar%f_annealing,error=error)
     IF (simpar%f_annealing /= 1.0_dp) simpar % annealing = .TRUE.
     CALL section_vals_val_get(md_section,"ENSEMBLE",i_val=simpar%ensemble,error=error)
     CALL section_vals_val_get(md_section,"STEPS",i_val=simpar%nsteps,error=error)
     CALL section_vals_val_get(md_section,"TIMESTEP",r_val=simpar%dt,error=error)
     CALL section_vals_val_get(md_section,"TEMPERATURE",r_val=simpar%temp_ext,error=error)
     CALL section_vals_val_get(md_section,"TEMP_TOL",r_val=simpar%temp_tol,error=error)
     CALL section_vals_val_get(md_section,"RESTART",i_val=simpar%read_type,error=error)
     CALL section_vals_val_get(md_section,"PV_AVA",l_val=simpar%pv_availability,error=error)
     CALL section_vals_val_get(md_section,"PV_NUM",l_val=simpar%pv_numerical,error=error)
     
     IF (simpar%ensemble == langevin_ensemble) THEN
       CALL section_vals_val_get(md_section, "GAMMA", r_val=simpar%gamma, error=error)
     END IF

     tmp_section => section_vals_get_subs_vals(globenv%input_file,"GLOBAL",error=error)
     CALL section_vals_val_get(globenv%input_file,"FORCE_EVAL%PRINT%TRAJECTORY%FORMAT",i_val=simpar%output_format,error=error)
     CALL section_vals_val_get(tmp_section,"WALLTIME",r_val=simpar%walltime,error=error)

     tmp_section => section_vals_get_subs_vals(globenv%input_file,"FORCE_EVAL%SUBSYS%CONSTRAINT",error=error)
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
!!     Initializes the simpar and mdio environment
!!
!!   AUTHOR
!!     teo
!!
!!   MODIFICATION HISTORY
!!  
!!
!!*** *************************************************************************  
  SUBROUTINE set_md_defaults(simpar, globenv, mdio)
    IMPLICIT NONE
    ! Arguments
    TYPE(simulation_parameters_type), &
         INTENT(OUT)                            :: simpar
    TYPE(global_environment_type), &
         INTENT(IN)                             :: globenv
    TYPE(mdio_parameters_type), INTENT(OUT)     :: mdio    
    ! Local Variables
    INTEGER                                  :: ia, ie

    simpar % units % time = 'FS'
    simpar % units % press = 'BAR'
    simpar % units % tau = 'FS'
    simpar % units % tau_cell = 'FS'
    simpar % units % verlet_skin = 'ANGSTROMS'
    simpar % output_format = dump_xmol
    simpar % ensemble = nve_ensemble
    simpar % nsteps = 0
    simpar % walltime = -1.0_dp ! default: don't care
    simpar % dt = 0.0_dp
    simpar % temp_ext = 0.0_dp
    simpar % temp_baro = 0.0_dp
    simpar % temp_tol = 0.0_dp
    simpar % p_ext = 0.0_dp
    simpar % f_annealing = 1.0_dp
    simpar % read_type = md_pos_restart
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
    simpar%md_method = Harris_MD_method
    simpar%mc_method = Harris_MC_method
    
    !..filenames
    CALL xstring ( globenv % project_name, ia, ie )
    mdio % crd_file_name   = globenv % project_name(ia:ie) // '.xyz'
    mdio % vel_file_name   = globenv % project_name(ia:ie) // '.vel'
    mdio % ptens_file_name = globenv % project_name(ia:ie) // '.ptens'
    mdio % ener_file_name  = globenv % project_name(ia:ie) // '.ener'
    mdio % pot_file_name   = globenv % project_name(ia:ie) // '.pot'
    mdio % temp_file_name  = globenv % project_name(ia:ie) // '.temp'
    mdio % dump_file_name  = globenv % project_name(ia:ie) // '.restart'
    mdio % coef_file_name  = globenv % project_name(ia:ie) // '.coef'
    mdio % hmc_file_name   = globenv % project_name(ia:ie) // '.hmc'

    !..print frequency
    mdio % icrd = 1
    mdio % icoef = 1
    mdio % ivel = 1
    mdio % iptens = 1
    mdio % iener = 1
    mdio % ipot = 1
    mdio % itemp = 1
    mdio % idump = 1
    mdio % iscreen = 1
    
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
!!       Type(global_environment_type), Intent (INOUT):: globenv
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

SUBROUTINE initialize_velocities ( simpar, part, globenv )


    TYPE(simulation_parameters_type), &
      INTENT(INOUT)                          :: simpar
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: part
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv

    INTEGER                                  :: i, iw, natoms
    REAL(KIND=dp)                                :: akin, denom, mass, sc, temp, mass_tot
    REAL(KIND=dp), DIMENSION(3)                  :: vcom, vcom_new
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(unit_convert_type), POINTER         :: unit

!------------------------------------------------------------------------------

  iw = globenv % scr
  IF (globenv % ionode .AND. globenv % print_level>1) THEN
     WRITE ( iw, '( A, A, A )' ) ' **********************', &
          ' begin of velocity initialization ', '***********************'
  END IF

! Initializing parameters
  vcom = 0.0_dp
  mass_tot = 0.0_dp
  vcom_new = 0.0_dp
  denom = 0.0_dp
  natoms = SIZE ( part )

  CALL get_cp2k_units ( unit ) 

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
          part ( i ) % v ( 1 ) = gasdev ( globenv % idum )
          part ( i ) % v ( 2 ) = gasdev ( globenv % idum )
          part ( i ) % v ( 3 ) = gasdev ( globenv % idum )
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

  IF ( (simpar % read_type == md_init_start ) .AND. &
        globenv % ionode .AND. globenv % print_level > 1 ) THEN
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

!******************************************************************************
