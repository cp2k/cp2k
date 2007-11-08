!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
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
  USE bibliography,                    ONLY: Evans1983,&
                                             Minary2003,&
                                             cite_reference
  USE cp2k_units,                      ONLY: get_cp2k_units,&
                                             unit_convert_type
  USE cp_linked_list_val,              ONLY: cp_sll_val_next,&
                                             cp_sll_val_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_generate_filename,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: &
       ensemble_name_types, isokin_ensemble, langevin_ensemble, nph_ensemble, &
       nph_uniaxial_damped_ensemble, nph_uniaxial_ensemble, npt_f_ensemble, &
       npt_i_ensemble, nve_ensemble
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_list_get,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE input_val_types,                 ONLY: val_get,&
                                             val_type
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE molecule_kind_types,             ONLY: fixd_constraint_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE parallel_rng_types,              ONLY: next_random_number
  USE particle_types,                  ONLY: particle_type
  USE shell_potential_types,           ONLY: shell_kind_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'md'

  ! Public data types

  PUBLIC :: simulation_parameters_type

  ! Public subroutines

  PUBLIC :: initialize_velocities,&
            read_md_section

  TYPE simulation_parameters_type
     INTEGER        :: nsteps
     REAL (KIND=dp) :: dt
     REAL (KIND=dp) :: temp_ext
     REAL (KIND=dp) :: temp_baro
     REAL (KIND=dp) :: temp_tol
     REAL (KIND=dp) :: comvel_tol
     REAL (KIND=dp) :: p_ext
     REAL (KIND=dp) :: cmass
     REAL (KIND=dp) :: v0
     REAL (KIND=dp) :: e0
     REAL (KIND=dp) :: v_shock
     REAL (KIND=dp) :: p0
     REAL (KIND=dp) :: f_annealing
     REAL (KIND=dp) :: gamma_nph
     INTEGER        :: ensemble
     LOGICAL :: constraint
     LOGICAL :: annealing
     LOGICAL :: dump_lm 
     INTEGER :: nfree, nfree_rot_transl
     INTEGER :: info_constraint
     INTEGER :: lagrange_multipliers
     REAL (KIND=dp) :: tau_cell
     REAL (KIND=dp) :: shake_tol
     ! Langevin Parameters
     REAL (KIND=dp) :: gamma
     REAL (KIND=dp) :: noisy_gamma
     REAL (KIND=dp) :: shadow_gamma
     REAL (KIND=dp) :: var_w
     ! RESPA Parameters
     LOGICAL        :: multi_time_switch, do_respa
     INTEGER        :: n_time_steps
     !SHELL parameters
     REAL (KIND=dp) :: temp_sh_ext
     REAL (KIND=dp) :: temp_sh_tol
     LOGICAL        :: temperature_per_kind
     LOGICAL        :: scale_temperature_per_kind
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
  SUBROUTINE read_md_section ( simpar,root_section, error)
    TYPE(simulation_parameters_type), &
      INTENT(OUT)                            :: simpar
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_md_section', &
      routineP = moduleN//':'//routineN

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

    motion_section => section_vals_get_subs_vals(root_section,"MOTION",error=error)
    md_section => section_vals_get_subs_vals(motion_section,"MD",error=error)

    logger => cp_error_get_logger(error)
    iw=cp_print_key_unit_nr(logger,md_section,"PRINT%PROGRAM_RUN_INFO",extension=".log",error=error)

    CALL read_md_low(simpar, root_section, md_section, error)
    IF ( iw > 0 ) WRITE (iw, *)

    ! Begin setup Langevin dynamics
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

    !..write some information to output
    IF ( iw > 0 ) THEN
       CALL get_cp2k_units ( unit )
       WRITE ( iw, '( A )' ) ' MD| Molecular Dynamics Protocol '
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
       IF (  simpar % ensemble == npt_i_ensemble .OR. simpar % ensemble == npt_f_ensemble ) THEN
          WRITE ( iw, '( A,A3,A, T71, F10.2 )' ) &
               ' MD| Pressure [','Bar','] ', simpar % p_ext*unit%bar
          WRITE ( iw, '( A,A4,A, T71, F10.2 )' ) &
               ' MD| Barostat time constant [','fs','] ', &
               simpar % tau_cell/unit%fs
       END IF
       IF (  simpar % ensemble == isokin_ensemble ) THEN
          CALL cite_reference(Evans1983)
          CALL cite_reference(Minary2003)
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
       ! Print warning for temp_tol
       IF (simpar % temp_tol > 0.0_dp) THEN
          CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
               "A temperature tolerance (TEMP_TOL) is used during the MD. "//&
               "Due to the velocity rescaling algorithm jumps may appear in the conserved quantity. "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
       END IF
    END IF

    CALL cp_print_key_finished_output(iw,logger,md_section,&
               "PRINT%PROGRAM_RUN_INFO", error=error)

  END SUBROUTINE read_md_section

!******************************************************************************
!!****** md/read_md_low [1.0] *
!!
!!   NAME
!!     read_md_low
!!
!!
!!   FUNCTION
!!     Low Level: Parses the MD input section
!!
!!   AUTHOR
!!     teo
!!
!!   MODIFICATION HISTORY
!!
!!
!!*** *************************************************************************
  SUBROUTINE read_md_low(simpar, root_section, md_section, error)
    TYPE(simulation_parameters_type), &
      INTENT(INOUT)                          :: simpar
    TYPE(section_vals_type), POINTER         :: root_section, md_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_md_low', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: explicit, failure
    TYPE(section_vals_type), POINTER         :: tmp_section

    failure = .FALSE.
    IF (.NOT.failure) THEN
       NULLIFY(tmp_section)
       CALL section_vals_val_get(md_section,"ANNEALING",r_val=simpar%f_annealing,error=error)
       IF (simpar%f_annealing /= 1.0_dp) simpar % annealing = .TRUE.
       CALL section_vals_val_get(md_section,"ENSEMBLE",i_val=simpar%ensemble,error=error)
       CALL section_vals_val_get(md_section,"STEPS",i_val=simpar%nsteps,error=error)
       CALL section_vals_val_get(md_section,"TIMESTEP",r_val=simpar%dt,error=error)
       CALL section_vals_val_get(md_section,"TEMPERATURE",r_val=simpar%temp_ext,error=error)
       CALL section_vals_val_get(md_section,"TEMP_TOL",r_val=simpar%temp_tol,error=error)
       CALL section_vals_val_get(md_section,"COMVEL_TOL",r_val=simpar%comvel_tol,error=error)
       CALL section_vals_val_get(md_section,"TEMP_KIND",l_val=simpar%temperature_per_kind,error=error)
       CALL section_vals_val_get(md_section,"SCALE_TEMP_KIND",l_val=simpar%scale_temperature_per_kind,error=error)

       IF (simpar%ensemble == langevin_ensemble) THEN
          CALL section_vals_val_get(md_section, "GAMMA", r_val=simpar%gamma, error=error)
          CALL section_vals_val_get(md_section, "NOISY_GAMMA", r_val=simpar%noisy_gamma, &
                                    error=error)
          CALL section_vals_val_get(md_section, "SHADOW_GAMMA", r_val=simpar%shadow_gamma, &
                                    error=error)
       END IF

       tmp_section => section_vals_get_subs_vals(root_section,"GLOBAL",error=error)
       tmp_section => section_vals_get_subs_vals(root_section,&
            "MOTION%CONSTRAINT",error=error)
       CALL section_vals_get(tmp_section,explicit=explicit,error=error)
       IF (explicit) THEN
          simpar%constraint = .TRUE.
          CALL section_vals_val_get(tmp_section,"SHAKE",r_val=simpar%shake_tol,error=error)
       END IF
       
       IF (simpar%ensemble == nph_uniaxial_ensemble.OR.simpar%ensemble == nph_uniaxial_damped_ensemble) THEN
          tmp_section => section_vals_get_subs_vals(md_section,"MSST",error=error)
          CALL section_vals_val_get(tmp_section,"PRESSURE",r_val=simpar%p0,error=error)
          CALL section_vals_val_get(tmp_section,"ENERGY",r_val=simpar%e0,error=error)
          CALL section_vals_val_get(tmp_section,"VOLUME",r_val=simpar%v0,error=error)
          CALL section_vals_val_get(tmp_section,"GAMMA",r_val=simpar%gamma_nph,error=error)
          IF ( simpar % gamma_nph /= 0.0_dp ) simpar % ensemble = nph_uniaxial_damped_ensemble
          CALL section_vals_val_get(tmp_section,"CMASS",r_val=simpar%cmass,error=error)
          CALL section_vals_val_get(tmp_section,"VSHOCK",r_val=simpar%v_shock,error=error)
       END IF
       
       SELECT CASE (simpar%ensemble)
       CASE(nph_ensemble, nph_uniaxial_damped_ensemble, nph_uniaxial_ensemble, npt_f_ensemble, npt_i_ensemble)
          tmp_section => section_vals_get_subs_vals(md_section,"BAROSTAT",error=error)
          CALL section_vals_val_get(tmp_section,"PRESSURE",r_val=simpar%p_ext,error=error)
          CALL section_vals_val_get(tmp_section,"TIMECON",r_val=simpar%tau_cell,error=error)
       END SELECT

       ! RESPA
       tmp_section => section_vals_get_subs_vals(md_section,"RESPA",error=error)
       CALL section_vals_get(tmp_section, explicit=simpar%do_respa, error=error)
       CALL section_vals_val_get(tmp_section,"FREQUENCY",i_val=simpar%n_time_steps,&
            ignore_required=.TRUE.,error=error)
       simpar%multi_time_switch = .FALSE.
       IF (simpar%do_respa) simpar%multi_time_switch = .TRUE.

       ! CORE-SHELL MODEL
       tmp_section => section_vals_get_subs_vals(md_section,"SHELL",error=error)
       CALL section_vals_get(tmp_section, explicit=explicit, error=error)
       IF (explicit) THEN
          CALL section_vals_val_get(tmp_section,"TEMPERATURE",r_val=simpar%temp_sh_ext,error=error)
          CALL section_vals_val_get(tmp_section,"TEMP_TOL",r_val=simpar%temp_sh_tol,error=error)
       END IF
    END IF

  END SUBROUTINE read_md_low

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

    simpar % ensemble = nve_ensemble
    simpar % nsteps = 0
    simpar % dt = 0.0_dp
    simpar % temp_ext = 0.0_dp
    simpar % temp_baro = 0.0_dp
    simpar % temp_tol = 0.0_dp
    simpar % comvel_tol = 0.0_dp
    simpar % p_ext = 0.0_dp
    simpar % f_annealing = 1.0_dp
    simpar % constraint = .FALSE.
    simpar % annealing = .FALSE.
    simpar % shake_tol = 1.0E-6_dp
    simpar % tau_cell = 1000.0_dp
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
    ! SHELL parameters
    simpar % temp_sh_ext = 0.0_dp
    simpar % temp_sh_tol = 0.0_dp
  END SUBROUTINE set_md_defaults

!******************************************************************************
!!****** md/compute_vcom [1.0] *
!!
!!   NAME
!!     compute_vcom
!!
!!   FUNCTION
!!     compute center of mass velocity
!!     *** is only used by initialize_velocities below ***
!!
!!   AUTHOR
!!     Toon Verstraelen <Toon.Verstraelen@gmail.com>
!!
!!   MODIFICATION HISTORY
!!     2007-11-6: created
!!
!!*** *************************************************************************
  SUBROUTINE compute_vcom(part,is_fixed,vcom,ecom)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    LOGICAL, DIMENSION(:), INTENT(IN)        :: is_fixed
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: vcom
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: ecom

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_vcom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: denom, mass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    vcom  = 0.0_dp
    denom = 0.0_dp
    DO i = 1, SIZE(part)
       atomic_kind => part(i)%atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
       IF (.NOT.is_fixed(i)) THEN
          vcom = vcom + part(i)%v(:) * mass
          denom = denom + mass
       END IF
    END DO
    vcom = vcom/denom
    IF (PRESENT(ecom)) THEN
       ecom = 0.5_dp*denom*SUM(vcom*vcom)
    END IF

  END SUBROUTINE compute_vcom

!******************************************************************************
!!****** md/clone_core_shell_vel [1.0] *
!!
!!   NAME
!!     clone_core_shell_vel
!!
!!   FUNCTION
!!     Copy atom velocities into core and shell velocities
!!     *** is only used by initialize_velocities below ***
!!
!!   AUTHOR
!!     Toon Verstraelen <Toon.Verstraelen@gmail.com>
!!
!!   MODIFICATION HISTORY
!!     2007-11-6: created
!!
!!*** *************************************************************************
  SUBROUTINE clone_core_shell_vel(part,shell_part,core_part)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part, shell_part, core_part

    CHARACTER(len=*), PARAMETER :: routineN = 'clone_core_shell_vel', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: is_shell
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    DO i = 1, SIZE(part)
       atomic_kind => part(i)%atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind, shell_active=is_shell)
       IF (is_shell) THEN
          shell_part( part(i)%shell_index )%v(:) =  part(i)%v(:)
          core_part( part(i)%shell_index )%v(:)  =  part(i)%v(:)
       END IF
    END DO

  END SUBROUTINE clone_core_shell_vel

!******************************************************************************
!!****** md/compute_ekin [1.0] *
!!
!!   NAME
!!     compute_ekin
!!
!!   FUNCTION
!!     Compute the kinetic energy.
!!     *** is only used by initialize_velocities below ***
!!
!!   AUTHOR
!!     Toon Verstraelen <Toon.Verstraelen@gmail.com>
!!
!!   MODIFICATION HISTORY
!!     2007-11-6: created
!!
!!*** *************************************************************************
  FUNCTION compute_ekin(part) RESULT(ekin)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    REAL(KIND=dp)                            :: ekin

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_ekin', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: mass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    ekin = 0.0_dp
    DO i = 1, SIZE(part)
       atomic_kind => part(i)%atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
       ekin = ekin + 0.5_dp * mass * SUM(part(i)%v(:) * part(i)%v(:))
    END DO

  END FUNCTION compute_ekin

!******************************************************************************
!!****** md/rescale_vel [1.0] *
!!
!!   NAME
!!     rescale_vel
!!
!!   FUNCTION
!!     Rescale the velocity to mimic the given external kinetic temperature.
!!     Optionally also rescale vcom.
!!     *** is only used by initialize_velocities below ***
!!
!!   AUTHOR
!!     Toon Verstraelen <Toon.Verstraelen@gmail.com>
!!
!!   MODIFICATION HISTORY
!!     2007-11-6: created
!!
!!*** *************************************************************************
  SUBROUTINE rescale_vel(part,simpar,ekin,vcom)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    REAL(KIND=dp), INTENT(INOUT)             :: ekin
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT), OPTIONAL                :: vcom

    CHARACTER(len=*), PARAMETER :: routineN = 'rescale_vel', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: factor
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    factor = simpar%temp_ext / (2.0_dp * ekin ) * REAL(simpar%nfree,KIND=dp)
    ! Note:
    ! this rescaling is still wrong, it should take the masses into account
    ! rescaling is generally not correct, so needs fixing
    ekin = ekin * factor
    factor = SQRT(factor)
    DO i = 1, SIZE(part)
       atomic_kind => part(i)%atomic_kind
       part(i)%v(:) = factor*part(i)%v(:)
    END DO
    IF (PRESENT(vcom)) THEN
       vcom = factor*vcom
    END IF

  END SUBROUTINE rescale_vel

!******************************************************************************
!!****** md/subtract_vcom [1.0] *
!!
!!   NAME
!!     subtract_vcom
!!
!!   FUNCTION
!!     subtract center of mass velocity
!!     *** is only used by initialize_velocities below ***
!!
!!   AUTHOR
!!     Toon Verstraelen <Toon.Verstraelen@gmail.com>
!!
!!   MODIFICATION HISTORY
!!     2007-11-6: created
!!
!!*** *************************************************************************
  SUBROUTINE subtract_vcom(part,is_fixed,vcom)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    LOGICAL, DIMENSION(:), INTENT(IN)        :: is_fixed
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: vcom

    CHARACTER(len=*), PARAMETER :: routineN = 'subtract_vcom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    DO i = 1, SIZE(part)
       IF (.NOT.is_fixed(i)) THEN
          part(i)%v(:) = part(i)%v(:) - vcom(:)
       END IF
    END DO
  END SUBROUTINE subtract_vcom

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
!!     CJM,MK,Toon Verstraelen <Toon.Verstraelen@gmail.com>
!!
!!   MODIFICATION HISTORY
!!     - is_fixed removed from particle_type
!!     - 2007-11-07: Cleanup (TV)
!!
!!*** *************************************************************************
  SUBROUTINE initialize_velocities(simpar,part,para_env,globenv,molecule_kinds,&
       print_section, subsys_section, shell_present, shell_part, core_part, &
       force_rescaling, error)

    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(section_vals_type), POINTER         :: print_section, subsys_section
    LOGICAL, INTENT(IN)                      :: shell_present
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: shell_part, core_part
    LOGICAL, INTENT(IN)                      :: force_rescaling
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_velocities', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ifixd, &
                                                imolecule_kind, iw, natoms, &
                                                nfixd, nshell, shell_index, &
                                                stat
    LOGICAL                                  :: atomvel_explicit, &
                                                corevel_explicit, failure, &
                                                is_ok, shellvel_explicit
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: is_fixed
    REAL(KIND=dp)                            :: ecom, ekin, fac_massc, &
                                                fac_masss, mass, mass_tot, &
                                                temp
    REAL(KIND=dp), DIMENSION(3)              :: v, vc, vcom, vs
    REAL(KIND=dp), DIMENSION(:), POINTER     :: vel
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_sll_val_type), POINTER           :: atom_list, core_list, &
                                                shell_list
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(section_vals_type), POINTER         :: atomvel_section, &
                                                corevel_section, &
                                                shellvel_section
    TYPE(shell_kind_type), POINTER           :: shell
    TYPE(unit_convert_type), POINTER         :: unit
    TYPE(val_type), POINTER                  :: val

    CALL timeset(routineN,"I","",handle)

    ! Initializing parameters
    failure = .FALSE.
    NULLIFY (atomic_kind, fixd_list, atom_list, shell_list, core_list, logger, molecule_kind)
    NULLIFY (molecule_kind_set, shell, unit, val)

    NULLIFY (atomvel_section,shellvel_section, corevel_section)
    atomvel_section  => section_vals_get_subs_vals(subsys_section,"VELOCITY",error=error)
    shellvel_section => section_vals_get_subs_vals(subsys_section,"SHELL_VELOCITY",error=error)
    corevel_section => section_vals_get_subs_vals(subsys_section,"CORE_VELOCITY",error=error)
    
    ! logging
    logger => cp_error_get_logger(error)
    iw=cp_print_key_unit_nr(logger,print_section,"PROGRAM_RUN_INFO",extension=".log",error=error)
    IF (iw>0) THEN
       WRITE ( iw, '( A, A, A )' ) ' **********************', &
            ' begin of velocity initialization ', '***********************'
    END IF

    ! Initializing parameters
    natoms = SIZE(part)
    nshell = 0

    ! Core-Shell Model
    IF (shell_present) THEN
       CPPostcondition(ASSOCIATED(core_part),cp_failure_level,routineP,error,failure)
       CPPostcondition(ASSOCIATED(shell_part),cp_failure_level,routineP,error,failure)
       nshell = SIZE(shell_part)
    END IF

    ! Build a list of all fixed atoms (if any)
    ALLOCATE (is_fixed(natoms),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    is_fixed(:) = .FALSE.
    molecule_kind_set => molecule_kinds%els
    DO imolecule_kind=1,molecule_kinds%n_els
       molecule_kind => molecule_kind_set(imolecule_kind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,nfixd=nfixd,fixd_list=fixd_list)
       DO ifixd=1,nfixd
          IF (.NOT.fixd_list(ifixd)%restraint%active) is_fixed(fixd_list(ifixd)%fixd) = .TRUE.
       END DO
    END DO
    
    ! Compute the total mass when needed
    IF ( simpar % ensemble == nph_ensemble .OR.&
         simpar % ensemble == nph_uniaxial_ensemble .OR.&
         simpar % ensemble == nph_uniaxial_damped_ensemble ) THEN
       mass_tot = 0.0_dp
       DO i = 1, natoms
          atomic_kind => part(i)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
          mass_tot = mass_tot + mass
       END DO
       simpar % v_shock = simpar % v_shock * SQRT ( mass_tot )
    END IF

    ! Read or initialize the particle velocities
    CALL get_cp2k_units ( unit )
    CALL section_vals_get(atomvel_section,explicit=atomvel_explicit,error=error)
    CALL section_vals_get(shellvel_section,explicit=shellvel_explicit,error=error)
    CALL section_vals_get(corevel_section,explicit=corevel_explicit,error=error)
    CPPostcondition(shellvel_explicit.EQV.corevel_explicit,cp_failure_level,routineP,error,failure)

    IF (atomvel_explicit) THEN
       ! Read the atom velocities if explicitely given in the input file
       CALL section_vals_list_get(atomvel_section,"_DEFAULT_KEYWORD_",list=atom_list,error=error)
       DO i = 1, natoms
         is_ok=cp_sll_val_next(atom_list,val,error=error)
         CALL val_get(val,r_vals=vel,error=error)
         part(i)%v = vel
         IF (is_fixed(i)) part(i)%v = 0.0_dp
       END DO
       IF(shell_present) THEN 
          IF (shellvel_explicit) THEN
             ! If the atoms positions are given (?) and core and shell velocities are
             ! present in the input, read the latter.
             CALL section_vals_list_get(shellvel_section,"_DEFAULT_KEYWORD_",list=shell_list,error=error)
             CALL section_vals_list_get(corevel_section,"_DEFAULT_KEYWORD_",list=core_list,error=error)
             DO i = 1 , nshell
                is_ok=cp_sll_val_next(shell_list,val,error=error)
                CALL val_get(val,r_vals=vel,error=error)
                shell_part(i)%v  = vel
                is_ok=cp_sll_val_next(core_list,val,error=error)
                CALL val_get(val,r_vals=vel,error=error)
                core_part(i)%v  = vel
             END DO
          ELSE
             ! Otherwise, just copy atom velocties into shell and core velocities.
             CALL clone_core_shell_vel(part,shell_part,core_part)
          END IF
       END IF
       
       ! compute vcom, ecom and ekin
       CALL compute_vcom(part,is_fixed,vcom,ecom)
       ekin = compute_ekin(part) - ecom
       
       IF (simpar % nfree /= 0 .AND. force_rescaling) THEN
          CALL rescale_vel(part,simpar,ekin,vcom)
          ! After rescaling, the core and shell velocities must also adapt.
          DO i = 1, natoms
            shell_index = part(i)%shell_index
            IF(shell_present .AND. shell_index/=0) THEN
                atomic_kind => part(i)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass, shell=shell)
                fac_masss = shell%mass_shell/mass
                fac_massc = shell%mass_core/mass
                vs = shell_part(shell_index)%v
                vc = core_part(shell_index)%v
                
                shell_part(shell_index)%v(1) = part(i)%v(1) + fac_massc*(vs(1)-vc(1))
                shell_part(shell_index)%v(2) = part(i)%v(2) + fac_massc*(vs(2)-vc(2))
                shell_part(shell_index)%v(3) = part(i)%v(3) + fac_massc*(vs(3)-vc(3))
                core_part(shell_index)%v(1)  = part(i)%v(1) + fac_masss*(vc(1)-vs(1))
                core_part(shell_index)%v(2)  = part(i)%v(2) + fac_masss*(vc(2)-vs(2))
                core_part(shell_index)%v(3)  = part(i)%v(3) + fac_masss*(vc(3)-vs(3))
            END IF
          END DO
       END IF
    ELSE
       ! Initializing velocities deterministically on all processors, if not given in input
       DO i = 1, natoms
          atomic_kind => part(i)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
          IF (is_fixed(i)) THEN
             part(i)%v(1) = 0.0_dp
             part(i)%v(2) = 0.0_dp
             part(i)%v(3) = 0.0_dp
          ELSE
             part(i)%v(1) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
             part(i)%v(2) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
             part(i)%v(3) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
          END IF
       END DO
       
       ! subtract the vcom and compute ekin
       CALL compute_vcom(part,is_fixed,vcom)
       CALL subtract_vcom(part,is_fixed,vcom)
       ekin = compute_ekin(part)

       ! rescale the velocities
       CALL rescale_vel(part,simpar,ekin)
       
       ! Initialize the core and the shell velocity. Atom velocities are just
       ! copied so that the initial relative core-shell velocity is zero.
       IF (shell_present) THEN
         CALL clone_core_shell_vel(part, shell_part, core_part)
       ENDIF 
    END IF
    
    ! Recompute vcom, ecom and ekin for IO
    CALL compute_vcom(part,is_fixed,vcom,ecom)
    ekin = compute_ekin(part)

    IF (iw>0) THEN
       IF (simpar % nfree == 0) THEN
          CPPostcondition(ekin==0.0_dp,cp_failure_level,routineP,error,failure)
          temp = 0.0_dp
       ELSE
          temp = 2.0_dp * ekin / REAL ( simpar % nfree,KIND=dp)
       END IF
       WRITE ( iw, '( A, T61, F18.2, A2 )' ) ' Initial Temperature ', temp / unit%kelvin, " K"
       WRITE ( iw, '( A, T61, F20.12 )' ) ' Centre of mass velocity in direction x:', vcom ( 1 )
       WRITE ( iw, '( A, T61, F20.12 )' ) ' Centre of mass velocity in direction y:', vcom ( 2 )
       WRITE ( iw, '( A, T61, F20.12 )' ) ' Centre of mass velocity in direction z:', vcom ( 3 )
       WRITE ( iw, '( A, A, A, / )' ) ' ***********************', ' end of velocity initialization ', '************************'
    END IF
    
    DEALLOCATE (is_fixed,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    CALL cp_print_key_finished_output(iw,logger,print_section,"PROGRAM_RUN_INFO", error=error)
    
    CALL timestop(0.0_dp,handle)
    
  END SUBROUTINE initialize_velocities
  ! ***************************************************************************
  
END MODULE md
