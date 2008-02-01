!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
MODULE md
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE bibliography,                    ONLY: Evans1983,&
                                             Minary2003,&
                                             cite_reference
  USE cell_types,                      ONLY: cell_type
  USE cp2k_units,                      ONLY: get_cp2k_units,&
                                             unit_convert_type
  USE cp_linked_list_val,              ONLY: cp_sll_val_next,&
                                             cp_sll_val_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_generate_filename,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: &
       ensemble_name_types, isokin_ensemble, langevin_ensemble, nph_ensemble, &
       nph_uniaxial_damped_ensemble, nph_uniaxial_ensemble, npt_f_ensemble, &
       npt_i_ensemble, nve_ensemble, use_perd_none, use_perd_x, use_perd_xy, &
       use_perd_xyz, use_perd_xz, use_perd_y, use_perd_yz, use_perd_z
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_list_get,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE input_val_types,                 ONLY: val_get,&
                                             val_type
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE mathlib,                         ONLY: diamat_all
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

! *****************************************************************************
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
     LOGICAL :: angvel_zero 
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

CONTAINS

! *****************************************************************************
!> \brief reads the md section and prints out what has been read
!> \par History
!>      none
!> \author Teo, CJM
! *****************************************************************************
  SUBROUTINE read_md_section ( simpar,root_section, error)
    TYPE(simulation_parameters_type), &
      INTENT(OUT)                            :: simpar
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_md_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER                                  :: iprint, iw
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
       CALL section_vals_val_get(print_key,"EACH%MD",i_val=iprint,error=error)
       WRITE ( iw, '( A,T63,i10,A )' ) ' MD| Print MD information every',iprint, ' step(s)'
       WRITE ( iw, '( A,T20,A,T71,A10 )' ) ' MD| File type','Print frequency[steps]', 'File names'

       print_key => section_vals_get_subs_vals(motion_section,"PRINT%TRAJECTORY",error=error)
       CALL section_vals_val_get(print_key,"EACH%MD",i_val=iprint,error=error)
       filename = cp_print_key_generate_filename(logger,print_key,&
            extension=".xyz",middle_name="pos",my_local=.FALSE.,error=error)
       WRITE ( iw, '( A,T20,i10,T31,A50 )' ) ' MD| Coordinates',iprint, &
            ADJUSTR ( TRIM(filename) )

       IF ( simpar % ensemble == nph_ensemble .OR.&
            simpar % ensemble == nph_uniaxial_ensemble .OR.&
            simpar % ensemble == nph_uniaxial_damped_ensemble.OR.&
            simpar % ensemble == npt_i_ensemble.OR.&
            simpar % ensemble == npt_f_ensemble ) THEN

          print_key => section_vals_get_subs_vals(motion_section,"PRINT%CELL",error=error)
          CALL section_vals_val_get(print_key,"EACH%MD",i_val=iprint,error=error)
          filename = cp_print_key_generate_filename(logger,print_key,&
               extension=".cell",my_local=.FALSE.,error=error)
          WRITE ( iw, '( A,T20,i10,T31,A50 )' ) ' MD| Simulation Cell',iprint, &
               ADJUSTR ( TRIM(filename) )
       END IF

       print_key => section_vals_get_subs_vals(motion_section,"PRINT%VELOCITIES",error=error)
       CALL section_vals_val_get(print_key,"EACH%MD",i_val=iprint,error=error)
       filename = cp_print_key_generate_filename(logger,print_key,&
            extension=".xyz",middle_name="vel",my_local=.FALSE.,error=error)
       WRITE ( iw, '( A,T20,i10,T31,A50 )' ) ' MD| Velocities',iprint,&
            ADJUSTR ( TRIM(filename) )

       print_key => section_vals_get_subs_vals(motion_section,"MD%PRINT%ENERGY",error=error)
       CALL section_vals_val_get(print_key,"EACH%MD",i_val=iprint,error=error)
       filename =cp_print_key_generate_filename(logger,print_key,&
            extension=".ener",my_local=.FALSE.,error=error)
       WRITE ( iw, '( A,T20,i10,T31,A50 )' ) ' MD| Energies',iprint, &
            ADJUSTR ( TRIM(filename) )

       print_key => section_vals_get_subs_vals(motion_section,"PRINT%RESTART",error=error)
       CALL section_vals_val_get(print_key,"EACH%MD",i_val=iprint,error=error)
       filename = cp_print_key_generate_filename(logger,print_key,&
            extension=".restart",my_local=.FALSE.,error=error)
       WRITE ( iw, '( A,T20,i10,T31,A50 )' ) ' MD| Dump',iprint, &
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

! *****************************************************************************
!> \brief Low Level: Parses the MD input section
!> \author teo
! *****************************************************************************
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
       CALL section_vals_val_get(md_section,"ANGVEL_ZERO",l_val=simpar%angvel_zero,error=error)
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

! *****************************************************************************
!> \brief Initializes the simpar
!> \author teo
! *****************************************************************************
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
    simpar % angvel_zero = .FALSE.
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

! *****************************************************************************
!> \brief compute center of mass position
!>      *** is only used by initialize_velocities below ***
!> \par History
!>      2007-11-6: created
!> \author Toon Verstraelen <Toon.Verstraelen@gmail.com>
! *****************************************************************************
  SUBROUTINE compute_rcom(part,is_fixed,rcom)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    INTEGER, DIMENSION(:), INTENT(IN)        :: is_fixed
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: rcom

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_rcom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: denom, mass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    rcom(:)  = 0.0_dp
    denom = 0.0_dp
    DO i = 1, SIZE(part)
       atomic_kind => part(i)%atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
       SELECT CASE(is_fixed(i))
       CASE(use_perd_x,use_perd_y,use_perd_z,use_perd_xy,use_perd_xz,use_perd_yz,use_perd_none)
          rcom(1) = rcom(1) + part(i)%r(1) * mass
          rcom(2) = rcom(2) + part(i)%r(2) * mass
          rcom(3) = rcom(3) + part(i)%r(3) * mass
          denom = denom + mass
       END SELECT
    END DO
    rcom = rcom/denom

  END SUBROUTINE compute_rcom

! *****************************************************************************
!> \brief compute center of mass velocity
!>      *** is only used by initialize_velocities below ***
!> \par History
!>      2007-11-6: created
!> \author Toon Verstraelen <Toon.Verstraelen@gmail.com>
! *****************************************************************************
  SUBROUTINE compute_vcom(part,is_fixed,vcom,ecom)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    INTEGER, DIMENSION(:), INTENT(IN)        :: is_fixed
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
       SELECT CASE(is_fixed(i))
       CASE(use_perd_x,use_perd_y,use_perd_z,use_perd_xy,use_perd_xz,use_perd_yz,use_perd_none)
          vcom(1) = vcom(1) + part(i)%v(1) * mass
          vcom(2) = vcom(2) + part(i)%v(2) * mass
          vcom(3) = vcom(3) + part(i)%v(3) * mass
          denom = denom + mass
       END SELECT
    END DO
    vcom = vcom/denom
    IF (PRESENT(ecom)) THEN
       ecom = 0.5_dp*denom*SUM(vcom*vcom)
    END IF

  END SUBROUTINE compute_vcom

! *****************************************************************************
!> \brief Copy atom velocities into core and shell velocities
!>      *** is only used by initialize_velocities below ***
!> \par History
!>      2007-11-6: created
!> \author Toon Verstraelen <Toon.Verstraelen@gmail.com>
! *****************************************************************************
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

! *****************************************************************************
!> \brief Compute the kinetic energy. Does not subtract the center of mass kinetic
!>      energy.
!>      *** is only used by initialize_velocities below ***
!> \par History
!>      2007-11-6: created
!> \author Toon Verstraelen <Toon.Verstraelen@gmail.com>
! *****************************************************************************
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

! *****************************************************************************
!> \brief Rescale the velocity to mimic the given external kinetic temperature.
!>      Optionally also rescale vcom.
!>      *** is only used by initialize_velocities below ***
!> \par History
!>      2007-11-6: created
!> \author Toon Verstraelen <Toon.Verstraelen@gmail.com>
! *****************************************************************************
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

    factor = simpar%temp_ext / ( 2.0_dp * ekin ) * REAL(simpar%nfree,KIND=dp)
    ! Note:
    ! this rescaling is still wrong, it should take the masses into account
    ! rescaling is generally not correct, so needs fixing
    ekin = ekin * factor
    factor = SQRT(factor)
    DO i = 1, SIZE(part)
       part(i)%v(:) = factor*part(i)%v(:)
    END DO
    IF (PRESENT(vcom)) THEN
       vcom = factor*vcom
    END IF

  END SUBROUTINE rescale_vel

! *****************************************************************************
!> \brief subtract center of mass velocity
!>      *** is only used by initialize_velocities below ***
!> \par History
!>      2007-11-6: created
!> \author Toon Verstraelen <Toon.Verstraelen@gmail.com>
! *****************************************************************************
  SUBROUTINE subtract_vcom(part,is_fixed,vcom)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    INTEGER, DIMENSION(:), INTENT(IN)        :: is_fixed
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: vcom

    CHARACTER(len=*), PARAMETER :: routineN = 'subtract_vcom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    DO i = 1, SIZE(part)
       SELECT CASE(is_fixed(i))
       CASE(use_perd_x)
          part(i)%v(2) = part(i)%v(2) - vcom(2)
          part(i)%v(3) = part(i)%v(3) - vcom(3)
       CASE(use_perd_y)
          part(i)%v(1) = part(i)%v(1) - vcom(1)
          part(i)%v(3) = part(i)%v(3) - vcom(3)
       CASE(use_perd_z)
          part(i)%v(1) = part(i)%v(1) - vcom(1)
          part(i)%v(2) = part(i)%v(2) - vcom(2)
       CASE(use_perd_xy)
          part(i)%v(3) = part(i)%v(3) - vcom(3)
       CASE(use_perd_xz)
          part(i)%v(2) = part(i)%v(2) - vcom(2)
       CASE(use_perd_yz)
          part(i)%v(1) = part(i)%v(1) - vcom(1)
       CASE(use_perd_none)
          part(i)%v(:) = part(i)%v(:) - vcom(:)
       END SELECT
    END DO
  END SUBROUTINE subtract_vcom

! *****************************************************************************
!> \brief compute the angular velocity
!>      *** is only used by initialize_velocities below ***
!> \par History
!>      2007-11-9: created
!> \author Toon Verstraelen <Toon.Verstraelen@gmail.com>
! *****************************************************************************
  SUBROUTINE compute_vang(part,is_fixed,rcom,vang)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    INTEGER, DIMENSION(:), INTENT(IN)        :: is_fixed
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rcom
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: vang

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_vang', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: mass, proj
    REAL(KIND=dp), DIMENSION(3)              :: evals, mang, r
    REAL(KIND=dp), DIMENSION(3, 3)           :: iner
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    mang(:)   = 0.0_dp
    iner(:,:) = 0.0_dp
    DO i=1,SIZE(part)
       ! compute angular momentum and inertia tensor
       SELECT CASE(is_fixed(i))
       CASE(use_perd_x,use_perd_y,use_perd_z,use_perd_xy,use_perd_xz,use_perd_yz,use_perd_none)
          r(:) = part(i)%r(:) - rcom(:)
          atomic_kind => part(i)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
          mang(1) = mang(1) + mass*(r(2)*part(i)%v(3) - r(3)*part(i)%v(2))
          mang(2) = mang(2) + mass*(r(3)*part(i)%v(1) - r(1)*part(i)%v(3))
          mang(3) = mang(3) + mass*(r(1)*part(i)%v(2) - r(2)*part(i)%v(1))

          iner(1,1) = iner(1,1) + mass*(r(2)*r(2) + r(3)*r(3))
          iner(2,2) = iner(2,2) + mass*(r(3)*r(3) + r(1)*r(1))
          iner(3,3) = iner(3,3) + mass*(r(1)*r(1) + r(2)*r(2))

          iner(1,2) = iner(1,2) - mass*r(1)*r(2)
          iner(2,3) = iner(2,3) - mass*r(2)*r(3)
          iner(3,1) = iner(3,1) - mass*r(3)*r(1)
       END SELECT
    END DO
    iner(2,1) = iner(1,2)
    iner(3,2) = iner(2,3)
    iner(1,3) = iner(3,1)

    ! Take the safest route, i.e. diagonalize the inertia tensor and solve
    ! the angular velocity only with the non-zero eigenvalues. A plain inversion
    ! would fail for linear molecules.
    CALL diamat_all(iner, evals)

    vang(:) = 0.0_dp
    DO i=1,3
       IF (evals(i) > 0.0_dp) THEN
         proj = SUM(iner(:,i)*mang)/evals(i)
         vang(1) = vang(1) + proj*iner(1,i)
         vang(2) = vang(2) + proj*iner(2,i)
         vang(3) = vang(3) + proj*iner(3,i)
       END IF
    END DO

  END SUBROUTINE compute_vang

! *****************************************************************************
!> \brief subtract the angular velocity
!>      *** is only used by initialize_velocities below ***
!> \par History
!>      2007-11-9: created
!> \author Toon Verstraelen <Toon.Verstraelen@gmail.com>
! *****************************************************************************
  SUBROUTINE subtract_vang(part,is_fixed,rcom,vang)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    INTEGER, DIMENSION(:), INTENT(IN)        :: is_fixed
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rcom, vang

    CHARACTER(len=*), PARAMETER :: routineN = 'subtract_vang', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(KIND=dp), DIMENSION(3)              :: r

    DO i=1,SIZE(part)
       r(:) = part(i)%r(:) - rcom(:)
       SELECT CASE(is_fixed(i))
       CASE(use_perd_x)
          part(i)%v(2) = part(i)%v(2) - (vang(3)*r(1) - vang(1)*r(3))
          part(i)%v(3) = part(i)%v(3) - (vang(1)*r(2) - vang(2)*r(1))
       CASE(use_perd_y)
          part(i)%v(1) = part(i)%v(1) - (vang(2)*r(3) - vang(3)*r(2))
          part(i)%v(3) = part(i)%v(3) - (vang(1)*r(2) - vang(2)*r(1))
       CASE(use_perd_z)
          part(i)%v(1) = part(i)%v(1) - (vang(2)*r(3) - vang(3)*r(2))
          part(i)%v(2) = part(i)%v(2) - (vang(3)*r(1) - vang(1)*r(3))
       CASE(use_perd_xy)
          part(i)%v(3) = part(i)%v(3) - (vang(1)*r(2) - vang(2)*r(1))
       CASE(use_perd_xz)
          part(i)%v(2) = part(i)%v(2) - (vang(3)*r(1) - vang(1)*r(3))
       CASE(use_perd_yz)
          part(i)%v(1) = part(i)%v(1) - (vang(2)*r(3) - vang(3)*r(2))
       CASE(use_perd_none)
          part(i)%v(1) = part(i)%v(1) - (vang(2)*r(3) - vang(3)*r(2))
          part(i)%v(2) = part(i)%v(2) - (vang(3)*r(1) - vang(1)*r(3))
          part(i)%v(3) = part(i)%v(3) - (vang(1)*r(2) - vang(2)*r(1))
       END SELECT
    END DO

  END SUBROUTINE subtract_vang

! *****************************************************************************
!> \brief Initializes the velocities to the Maxwellian distribution
!> \par History
!>      - is_fixed removed from particle_type
!>      - 2007-11-07: Cleanup (TV)
!>      - 2007-11-09: Added angvel_zero feature
!> \author CJM,MK,Toon Verstraelen <Toon.Verstraelen@gmail.com>
! *****************************************************************************
  SUBROUTINE initialize_velocities(simpar, part, force_env, para_env, globenv, &
       molecule_kinds, print_section, subsys_section, shell_present, shell_part, &
       core_part, force_rescaling, error)

    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    TYPE(force_env_type), POINTER            :: force_env
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
                                                nshell, shell_index, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: is_fixed
    LOGICAL                                  :: atomvel_explicit, &
                                                corevel_explicit, failure, &
                                                is_ok, shellvel_explicit
    REAL(KIND=dp)                            :: ecom, ekin, fac_massc, &
                                                fac_masss, mass, mass_tot, &
                                                temp
    REAL(KIND=dp), DIMENSION(3)              :: rcom, vang, vc, vcom, vs
    REAL(KIND=dp), DIMENSION(:), POINTER     :: vel
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
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
    
    ! Logging
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

    is_fixed = use_perd_none
    molecule_kind_set => molecule_kinds%els
    DO imolecule_kind=1,molecule_kinds%n_els
       molecule_kind => molecule_kind_set(imolecule_kind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,fixd_list=fixd_list)
       IF (ASSOCIATED(fixd_list)) THEN
          DO ifixd=1,SIZE(fixd_list)
             IF (.NOT.fixd_list(ifixd)%restraint%active) is_fixed(fixd_list(ifixd)%fixd) = fixd_list(ifixd)%itype
          END DO
       END IF
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
         SELECT CASE(is_fixed(i))
         CASE (use_perd_x)
            part(i)%v(1) = 0.0_dp
         CASE (use_perd_y)
            part(i)%v(2) = 0.0_dp
         CASE (use_perd_z)
            part(i)%v(3) = 0.0_dp
         CASE (use_perd_xy)
            part(i)%v(1) = 0.0_dp
            part(i)%v(2) = 0.0_dp
         CASE (use_perd_xz)
            part(i)%v(1) = 0.0_dp
            part(i)%v(3) = 0.0_dp
         CASE (use_perd_yz)
            part(i)%v(2) = 0.0_dp
            part(i)%v(3) = 0.0_dp
         CASE (use_perd_xyz)
            part(i)%v    = 0.0_dp
         END SELECT
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
          part(i)%v(1) = 0.0_dp
          part(i)%v(2) = 0.0_dp
          part(i)%v(3) = 0.0_dp
          SELECT CASE(is_fixed(i))
          CASE (use_perd_x)
             part(i)%v(2) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
             part(i)%v(3) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
          CASE (use_perd_y)
             part(i)%v(1) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
             part(i)%v(3) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
          CASE (use_perd_z)
             part(i)%v(1) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
             part(i)%v(2) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
          CASE (use_perd_xy)
             part(i)%v(3) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
          CASE (use_perd_xz)
             part(i)%v(2) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
          CASE (use_perd_yz)
             part(i)%v(1) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
          CASE (use_perd_none)
             part(i)%v(1) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
             part(i)%v(2) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
             part(i)%v(3) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
          END SELECT
       END DO
       
       ! Subtract the vcom
       CALL compute_vcom(part,is_fixed,vcom)
       CALL subtract_vcom(part,is_fixed,vcom)
       
       ! If requested and the system is not periodic, subtract the angular velocity
       CALL force_env_get(force_env, cell=cell, error=error)
       IF (SUM(cell%perd(1:3)) == 0 .AND. simpar%angvel_zero) THEN
          CALL compute_rcom(part,is_fixed,rcom)
          CALL compute_vang(part,is_fixed,rcom,vang)
          CALL subtract_vang(part,is_fixed,rcom,vang)
       END IF

       ! Rescale the velocities
       ekin = compute_ekin(part)
       CALL rescale_vel(part,simpar,ekin)
       
       ! Initialize the core and the shell velocity. Atom velocities are just
       ! copied so that the initial relative core-shell velocity is zero.
       IF (shell_present) THEN
         CALL clone_core_shell_vel(part, shell_part, core_part)
       ENDIF 
    END IF

    IF (iw>0) THEN
       ! Recompute vcom, ecom and ekin for IO
       CALL compute_vcom(part,is_fixed,vcom,ecom)
       ekin = compute_ekin(part) - ecom
       IF (simpar % nfree == 0) THEN
          CPPostcondition(ekin==0.0_dp,cp_failure_level,routineP,error,failure)
          temp = 0.0_dp
       ELSE
          temp = 2.0_dp * ekin / REAL ( simpar % nfree,KIND=dp)
       END IF
       WRITE (iw, '( A, T61, F18.2, A2 )' ) ' Initial Temperature ', temp / unit%kelvin, " K"
       WRITE (iw, '( A, T21, F20.12 , F20.12 , F20.12 )' ) ' COM velocity:', vcom ( 1 ), vcom ( 2 ), vcom ( 3 )

       ! compute and log rcom and vang if not periodic
       CALL force_env_get(force_env, cell=cell, error=error)
       IF (SUM(cell%perd(1:3)) == 0) THEN
          CALL compute_rcom(part,is_fixed,rcom)
          CALL compute_vang(part,is_fixed,rcom,vang)
          WRITE (iw, '( A, T21, F20.12 , F20.12 , F20.12 )' ) ' COM position:', rcom ( 1 ), rcom ( 2 ), rcom ( 3 )
          WRITE (iw, '( A, T21, F20.12 , F20.12 , F20.12 )' ) ' Angular velocity:', vang ( 1 ), vang ( 2 ), vang ( 3 )
       END IF
       WRITE ( iw, '( A, A, A, / )' ) ' ***********************', ' end of velocity initialization ', '************************'
    END IF
    
    DEALLOCATE (is_fixed,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    CALL cp_print_key_finished_output(iw,logger,print_section,"PROGRAM_RUN_INFO", error=error)
    CALL timestop(0.0_dp,handle)
    
  END SUBROUTINE initialize_velocities
  
END MODULE md
