!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/md_energies [1.0] *
!!
!!   NAME
!!     md_energies
!!
!!   FUNCTION
!!     prints all energy info per timestep to the screen or to
!!     user defined output files
!!
!!   AUTHOR
!!     Joost VandeVondele (copy from md_fist_energies)
!!
!!   MODIFICATIONS
!!     - New MD data are appended to the old data (15.09.2003,MK)
!!
!******************************************************************************

MODULE md_energies

  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell
  USE cp2k_units,                      ONLY: get_cp2k_units,&
                                             unit_convert_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type,&
                                             use_mixed_force
  USE input_constants,                 ONLY: dump_dcd,&
                                             dump_xmol,&
                                             ensemble_name_types,&
                                             nph_ensemble,&
                                             nph_uniaxial_damped_ensemble,&
                                             nph_uniaxial_ensemble,&
                                             npt_f_ensemble,&
                                             npt_i_ensemble,&
                                             reftraj_ensemble
  USE input_section_types,             ONLY: section_get_ival,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             sp
  USE machine,                         ONLY: m_flush
  USE md,                              ONLY: simulation_parameters_type
  USE md_conserved_quantities,         ONLY: compute_conserved_quantity
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_ener_type,&
                                             md_environment_type,&
                                             set_md_env,&
                                             zero_md_ener
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: write_particle_coordinates,&
                                             write_structure_data
  USE physcon,                         ONLY: angstrom,&
                                             femtoseconds,&
                                             kelvin
  USE qmmm_types,                      ONLY: qmmm_env_qm_type
  USE reftraj_types,                   ONLY: reftraj_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
  USE f77_blas
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'md_energies'

  TYPE average_quantities_type
     REAL ( dp ) :: avetemp, avepot, avekin, avevol, aveca, avecb, avecc
     REAL ( dp ) :: avetemp_baro, avehugoniot, avecpu
     REAL ( dp ) :: aveal, avebe, avega, avepress, avekinc, avetempc, avepxx
     REAL ( dp ) :: avetemp_qm, avekin_qm
  END TYPE average_quantities_type


  PUBLIC :: initialize_md_ener,&
            md_energy,&
            md_ener_reftraj,&
            md_write_output,&
            get_output_format

!!***
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** md_energies/initialize_md_ener [1.0] *
!!
!!   NAME
!!     initialize_md_ener
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     -10-2007 created 
!!
!!*****************************************************************************

  SUBROUTINE initialize_md_ener(md_ener,force_env,simpar,error)

    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'initialize_md_ener', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat, nkind
    LOGICAL                                  :: failure, shell_adiabatic
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles, shell_particles
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env

    failure =.FALSE.

    NULLIFY(qmmm_env, subsys)
    NULLIFY(atomic_kinds, atomic_kind_set, qmmm_env, particles, shell_particles)

    CPPrecondition(ASSOCIATED(md_ener),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)

    CALL force_env_get(force_env, subsys=subsys, qmmm_env=qmmm_env, error=error)
    CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds, particles=particles,&
         shell_particles=shell_particles,&
         error=error)
    atomic_kind_set => atomic_kinds%els
    nkind =  SIZE(atomic_kind_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               shell_adiabatic=shell_adiabatic)



    md_ener%nfree       = simpar%nfree
    md_ener%nfree_qm    = -HUGE(0)
    md_ener%nfree_shell = -HUGE(0)

    IF(shell_adiabatic) THEN
      md_ener%nfree_shell = 3*(shell_particles%n_els)
    END IF
    IF(ASSOCIATED(qmmm_env)) THEN
      ! The degrees of freedom for the quantum part of the system
      ! are set to 3*Number of QM atoms and to simpar%nfree in case all the MM
      ! system is treated at QM level (not really QM/MM, just for consistency).
      ! The degree of freedom will not be correct if 1-3 atoms are treated only
      ! MM. In this case we should take care of rotations
      md_ener%nfree_qm = 3*SIZE(qmmm_env%qm_atom_index)
      IF (md_ener%nfree_qm == 3*(particles%n_els)) md_ener%nfree_qm = md_ener%nfree
    END IF
    
    IF(simpar%temperature_per_kind) THEN
      ALLOCATE(md_ener%temp_kind(nkind), STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(md_ener%ekin_kind(nkind), STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(md_ener%nfree_kind(nkind),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      md_ener%nfree_kind = 0
      
      IF(shell_adiabatic) THEN 
        ALLOCATE(md_ener%temp_shell_kind(nkind), STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(md_ener%ekin_shell_kind(nkind), STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(md_ener%nfree_shell_kind(nkind), STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        md_ener%nfree_shell_kind = 0
      END IF

    END IF   
    CALL zero_md_ener(md_ener, tkind=simpar%temperature_per_kind, &
         tshell=shell_adiabatic, error=error)
    md_ener%epot = 0.0_dp
 

  END SUBROUTINE initialize_md_ener

!******************************************************************************
!!****** md_energies/md_energy [1.0] *
!!
!!   NAME
!!     md_energy
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     -10-2007 created 
!!
!!*****************************************************************************

  SUBROUTINE md_energy(md_env, md_ener, error)

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'md_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: natom
    LOGICAL                                  :: init, shell_adiabatic, &
                                                tcoeff, tkind, tshell
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell_ref
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar

    NULLIFY(atomic_kinds, atomic_kind_set, cell_ref, dyn_coeff_set, force_env,&
            particles, subsys, simpar)
    CALL get_md_env(md_env=md_env, init=init, force_env=force_env, &
         simpar=simpar, error=error)

    CALL force_env_get(force_env, cell_ref=cell_ref, dyn_coeff_set=dyn_coeff_set, &
         potential_energy=md_ener%epot, subsys=subsys,  error=error)
    
    CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds, error=error)
    atomic_kind_set => atomic_kinds%els
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               shell_adiabatic=shell_adiabatic)

    tcoeff = ASSOCIATED(dyn_coeff_set)
    tkind = simpar%temperature_per_kind
    tshell = shell_adiabatic

    CALL cp_subsys_get(subsys,particles=particles,error=error)
    natom=particles%n_els

    !---------------------------------------------
    ! Initialize V0 properly for nph calculations
    !---------------------------------------------
    IF (init) THEN
       IF ( simpar % ensemble == nph_ensemble .OR.&
            simpar % ensemble == nph_uniaxial_ensemble .OR.&
            simpar % ensemble == nph_uniaxial_damped_ensemble ) THEN
          IF ( simpar % v0 == 0.0_dp ) simpar % v0 = cell_ref % deth
       ENDIF
    ENDIF

    
    CALL compute_conserved_quantity(md_env, md_ener, tkind=tkind,&
         tcoeff=tcoeff, tshell=tshell, natom=natom, error=error)
    
  END SUBROUTINE md_energy

!!****f* md_energies/md_ener_reftraj [1.0] *
!!
!!   NAME
!!     md_ener_reftraj
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     10.2007 created 
!!
!!*** **********************************************************************

SUBROUTINE md_ener_reftraj(md_env,md_ener,error)
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'md_ener_reftraj', &
      routineP = moduleN//':'//routineN

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(reftraj_type), POINTER              :: reftraj


    CALL zero_md_ener(md_ener, tkind=.FALSE., tshell=.FALSE., error=error)

    CALL get_md_env(md_env=md_env,  force_env=force_env, reftraj=reftraj, error=error)
    
    IF(reftraj%info%eval_ef) THEN
      CALL force_env_get(force_env, potential_energy=md_ener%epot, error=error)
    ELSE
      md_ener%epot = reftraj%epot
      md_ener%delta_epot = (reftraj%epot - reftraj%epot0)/&
          REAL(reftraj%natom, kind=dp)*kelvin
    END IF
    
END SUBROUTINE  md_ener_reftraj

!******************************************************************************
!!****** md_energies/md_write_output [1.0] *
!!
!!   NAME
!!     md_write_output
!!
!!   FUNCTION
!!     This routine computes the conserved quantity, temperature
!!     and things like that and prints them out
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     - New MD data are appended to the old data (15.09.2003,MK)
!!
!!*****************************************************************************

  SUBROUTINE md_write_output(md_env, error)

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'md_write_output', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: fmd, my_act, my_ext, my_form, &
                                                my_pos, remark1, remark2, &
                                                title
    INTEGER :: box, cocrd, cofrc, covel, crd, delta_t, ene, ener_mix, frc, &
      handle, i, iskip, iw, nat, nkind, nshell, outformat, shcrd, shene, &
      shfrc, shvel, str, tempkind, tempshellkind, trsl, vel
    INTEGER, POINTER                         :: itimes
    LOGICAL :: failure, init, init_coforce, init_cotraj, init_covel, &
      init_force, init_shforce, init_shtraj, init_shvel, init_traj, init_vel, &
      is_mixed, qmmm, shell_adiabatic, shell_present
    REAL(dp)                                 :: abc( 3 ), cell_angle( 3 ), &
                                                econs, pv_scalar, pv_xx
    REAL(KIND=dp)                            :: harm_shell, hugoniot
    REAL(KIND=dp), POINTER                   :: constant, time, used_time
    REAL(KIND=dp), SAVE                      :: econs_avg, itimes_old, t_old
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(average_quantities_type), SAVE      :: averages
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE (reftraj_type), POINTER             :: reftraj
    TYPE(section_vals_type), POINTER         :: motion_section, root_section
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(unit_convert_type), POINTER         :: unit
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.
    NULLIFY(logger,unit)
    logger => cp_error_get_logger(error)
    CALL timeset(routineN,handle)
    CALL get_cp2k_units ( unit )

    ! Zeroing
    hugoniot = 0._dp
    econs = 0.0_dp
    shell_adiabatic = .FALSE.
    shell_present   = .FALSE.
    NULLIFY(dyn_coeff_set,motion_section, atomic_kinds, atomic_kind_set)
    NULLIFY(cell, subsys, force_env, md_ener, qmmm_env, reftraj)
    NULLIFY(core_particles, particles,shell_particles)
    NULLIFY(root_section,simpar,virial)

    CALL get_md_env(md_env=md_env, itimes=itimes, constant=constant,&
         t=time, used_time = used_time, simpar=simpar, force_env=force_env, &
         init=init, md_ener=md_ener, reftraj=reftraj, error=error)

    root_section   => force_env%root_section
    motion_section => section_vals_get_subs_vals(root_section,"MOTION",error=error)

    CALL force_env_get(force_env, cell=cell, subsys=subsys,&
         dyn_coeff_set=dyn_coeff_set, qmmm_env=qmmm_env, virial=virial, &
         error=error)
    qmmm = .FALSE.
    IF (md_ener%nfree_qm>0) qmmm = .TRUE.

    is_mixed= (force_env%in_use == use_mixed_force)

    CALL cp_subsys_get(subsys,particles=particles,error=error)
    nat=particles%n_els

    ! Computing the scalar pressure
    IF ( virial % pv_availability ) THEN
       pv_scalar = 0._dp
       DO i = 1, 3
          pv_scalar = pv_scalar + virial%pv_total(i,i)
       END DO
       pv_scalar = pv_scalar/3._dp/cell%deth
       pv_scalar = pv_scalar * unit % bar
       pv_xx = virial%pv_total(1,1)/cell%deth*unit%bar
    ENDIF

    CALL force_env_get(force_env,  qmmm_env=qmmm_env, error=error)

    CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds, error=error)
    atomic_kind_set => atomic_kinds%els
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                             shell_present=shell_present, &
                             shell_adiabatic=shell_adiabatic)

    CALL get_cell ( cell, abc = abc,  &
         alpha = cell_angle ( 3 ), &
         beta  = cell_angle ( 2 ), &
         gamma = cell_angle ( 1 ) )

    my_pos = "APPEND"
    my_act = "WRITE"
    IF (init.AND.(itimes==0)) THEN
       my_pos = "REWIND"
       my_act = "WRITE"
    END IF
    IF (init) THEN
       !-------------------------------------------------------------
       ! Performing protocol relevant to the first step of an MD run
       !-------------------------------------------------------------
       ! Print to the screen no matter what (e.g. ignore_should_output = .TRUE.)
       iw = cp_print_key_unit_nr(logger,motion_section,"MD%PRINT%PROGRAM_RUN_INFO",&
           extension=".mdLog",ignore_should_output=.TRUE.,file_position=my_pos,    &
           file_action=my_act,error=error)
       ! Computing the Hugoniot for NPH calculations
       IF ( simpar % ensemble == nph_ensemble .OR.&
            simpar % ensemble == nph_uniaxial_ensemble .OR.&
            simpar % ensemble == nph_uniaxial_damped_ensemble ) THEN
          IF ( simpar % e0 == 0._dp ) simpar % e0 = md_ener%epot + md_ener%ekin
          hugoniot = md_ener%epot + md_ener%ekin - simpar % e0 -  &
               .5_dp * ( pv_xx/unit%bar + simpar % p0 ) * &
               ( simpar % v0 - cell % deth )
       ENDIF
       ! Assigning the initial time-step and time
       itimes_old = itimes
       t_old      = time

       IF (simpar%ensemble==reftraj_ensemble) THEN
          itimes_old = reftraj%itimes0
          t_old = reftraj%time0
          reftraj%init = init
       END IF
       
       ! Write initial values of quantities of interest
       IF (iw>0) THEN
          WRITE (iw,*)
          WRITE (iw, '( A )' ) ' MD_ENERGIES| Initialization proceeding'
          WRITE (iw,*)
          WRITE ( iw, '( )' )
          WRITE ( iw, '( A,A )' ) ' ******************************** ', &
               'GO CP2K GO! **********************************'
          WRITE ( iw, '( A,A,T40,A,T60,1(1X,E20.12) )' ) ' INITIAL POTENTIAL ENERGY', &
               ADJUSTL( unit % e_label), '= ', md_ener%epot
          IF (simpar%ensemble/=reftraj_ensemble) THEN
            IF (.NOT.qmmm) THEN
              ! NO QM/MM info
              WRITE ( iw, '( A,A,T40,A,T60,1(1X,E20.12) )' ) ' INITIAL KINETIC ENERGY', &
                  ADJUSTL( unit % e_label), '= ', md_ener%ekin
              WRITE ( iw, '( A,A,T40,A,T60,1(1X,F20.3) )' ) ' INITIAL TEMPERATURE', &
                  ADJUSTL( unit % temp_label ), '= ', md_ener%temp_part
            ELSE
              WRITE ( iw, '( A,A,T40,A,T60,1(1X,E20.12) )' ) ' MM INITIAL KINETIC ENERGY', &
                  ADJUSTL( unit % e_label), '= ', md_ener%ekin
              WRITE ( iw, '( A,A,T40,A,T60,1(1X,E20.12) )' ) ' QM INITIAL KINETIC ENERGY', &
                  ADJUSTL( unit % e_label), '= ', md_ener%ekin_qm
              WRITE ( iw, '( A,A,T40,A,T60,1(1X,F20.3) )' ) ' MM INITIAL TEMPERATURE', &
                  ADJUSTL( unit % temp_label ), '= ', md_ener%temp_part
              WRITE ( iw, '( A,A,T40,A,T60,1(1X,F20.3) )' ) ' QM INITIAL TEMPERATURE', &
                  ADJUSTL( unit % temp_label ), '= ', md_ener%temp_qm
            END IF
          END IF
          IF ( simpar % ensemble == nph_ensemble .OR.&
               simpar % ensemble == nph_uniaxial_ensemble .OR.&
               simpar % ensemble == nph_uniaxial_damped_ensemble.OR.&
               simpar % ensemble == npt_i_ensemble.OR.&
               simpar % ensemble == npt_f_ensemble ) &
               WRITE ( iw, '( A,A,T40,A,T60,1(1X,E20.12) )' ) ' INITIAL BAROSTAT TEMP', &
               ADJUSTL( unit % temp_label), '= ', md_ener%temp_baro
          IF ( virial % pv_availability ) &
               WRITE ( iw, '( A,A,T40,A,T60,1(1X,E20.12) )' ) ' INITIAL PRESSURE', &
               ADJUSTL( unit % pres_label), '= ', pv_scalar
          IF ( simpar % ensemble == nph_ensemble .OR.&
               simpar % ensemble == nph_uniaxial_ensemble .OR.&
               simpar % ensemble == nph_uniaxial_damped_ensemble)&
               WRITE ( iw, '( A,A,T40,A,T60,1(1X,E20.12) )' ) ' INITIAL HUGONIOT CONSTRAINT', &
               ADJUSTL( unit % temp_label), '= ', hugoniot
          IF ( simpar % ensemble == nph_ensemble .OR.&
               simpar % ensemble == nph_uniaxial_ensemble .OR.&
               simpar % ensemble == nph_uniaxial_damped_ensemble)&
               WRITE ( iw, '( A,A,T40,A,T60,1(1X,E20.12) )' ) ' INITIAL E0', &
               ADJUSTL( unit % e_label), '= ', simpar % e0
          WRITE ( iw, '( A,A,T40,A,T60,1(1X,E20.12) )' ) ' INITIAL VOLUME', &
               ADJUSTL( unit % vol_label), '= ', cell%deth
          WRITE ( iw, '( A,A,T29,A,T33,3(1X,E15.7) )' ) ' INITIAL CELL LNTHS', &
               ADJUSTL(unit % l_label), '= ', abc(1), abc(2), abc(3)
          WRITE ( iw, '( A,A,T29,A,T33,3(1X,E15.7) )' ) ' INITIAL CELL ANGLS', &
               ADJUSTL(unit % angl_label), '= ', cell_angle(3), cell_angle(2), cell_angle(1)
          WRITE ( iw, '( A,A )' ) ' ******************************** ', &
               'GO CP2K GO! **********************************'
       ENDIF
       CALL cp_print_key_finished_output(iw,logger,motion_section,&
            "MD%PRINT%PROGRAM_RUN_INFO",ignore_should_output=.TRUE., error=error)
    ELSE
       !-------------------------------------------------------------
       ! Performing protocol for anything beyond the first step of MD
       !-------------------------------------------------------------
       IF ( simpar % ensemble == nph_ensemble .OR.&
            simpar % ensemble == nph_uniaxial_ensemble .OR.&
            simpar % ensemble == nph_uniaxial_damped_ensemble)&
            hugoniot = md_ener%epot + md_ener%ekin - simpar % e0 -  &
            .5_dp * ( pv_xx/unit%bar + simpar % p0 ) * &
            ( simpar % v0 - cell % deth )

       IF (simpar%ensemble==reftraj_ensemble) THEN
          time = reftraj%time
          econs = md_ener%delta_epot
          itimes = reftraj%itimes
       ELSE
          time = t_old + REAL(itimes-itimes_old,dp)*simpar%dt
          econs = md_ener%delta_cons!(md_ener%constant-constant)/nat*kelvin
       END IF

       ! average quantities:
       delta_t = itimes - itimes_old
       SELECT CASE (delta_t)
          CASE (1)
          IF ( virial % pv_availability ) averages%avepress = pv_scalar
          IF ( virial % pv_availability ) averages%avepxx = pv_xx
          averages%avecpu = used_time
          averages%avehugoniot = hugoniot
          averages%avetemp_baro = md_ener%temp_baro
          averages%avepot = md_ener%epot
          averages%avekin = md_ener%ekin
          averages%avetemp = md_ener%temp_part
          averages%avekin_qm  = md_ener%ekin_qm
          averages%avetemp_qm = md_ener%temp_qm
          averages%avevol = cell % deth
          averages%aveca = abc ( 1 )
          averages%avecb = abc ( 2 )
          averages%avecc = abc ( 3 )
          averages%aveal = cell_angle ( 3 )
          averages%avebe = cell_angle ( 2 )
          averages%avega = cell_angle ( 1 )
          econs_avg =  0._dp
       CASE DEFAULT
          IF ( virial % pv_availability ) &
               averages%avepress = (averages%avepress*REAL(delta_t-1,dp)+pv_scalar)/ &
               REAL ( delta_t, dp )
          IF ( virial % pv_availability ) &
               averages%avepxx = (averages%avepxx*REAL(delta_t-1,dp)+pv_xx)/ &
               REAL ( delta_t, dp )
          averages%avecpu = (averages%avecpu*REAL(delta_t-1,dp)+used_time)/ &
               REAL(delta_t,dp)
          averages%avehugoniot = (averages%avehugoniot*REAL(delta_t-1,dp)+hugoniot)/ &
               REAL(delta_t,dp)
          averages%avetemp_baro= (averages%avetemp_baro*REAL(delta_t-1,dp)+md_ener%temp_baro)/ &
               REAL(delta_t,dp)
          averages%avepot = (averages%avepot*REAL(delta_t-1,dp)+md_ener%epot)/ &
               REAL(delta_t,dp)
          averages%avekin = (averages%avekin*REAL(delta_t-1,dp)+md_ener%ekin)/ &
               REAL(delta_t,dp)
          averages%avetemp = (averages%avetemp*REAL(delta_t-1,dp)+md_ener%temp_part)/ &
               REAL(delta_t,dp)
          averages%avekin_qm = (averages%avekin_qm*REAL(delta_t-1,dp)+md_ener%ekin_qm)/ &
               REAL(delta_t,dp)
          averages%avetemp_qm = (averages%avetemp_qm*REAL(delta_t-1,dp)+md_ener%temp_qm)/ &
               REAL(delta_t,dp)
          averages%avevol = (averages%avevol*REAL(delta_t-1,dp)+cell%deth)/ &
               REAL(delta_t,dp)
          averages%aveca = (averages%aveca*REAL(delta_t-1,dp)+abc(1))/REAL(delta_t,dp)
          averages%avecb = (averages%avecb*REAL(delta_t-1,dp)+abc(2))/REAL(delta_t,dp)
          averages%avecc = (averages%avecc*REAL(delta_t-1,dp)+abc(3))/REAL(delta_t,dp)
          averages%aveal = (averages%aveal*REAL(delta_t-1,dp)+cell_angle(3))/ &
               REAL(delta_t,dp)
          averages%avebe = (averages%avebe*REAL(delta_t-1,dp)+cell_angle(2))/ &
               REAL(delta_t,dp)
          averages%avega = (averages%avega*REAL(delta_t-1,dp)+cell_angle(1))/ &
               REAL(delta_t,dp)
          econs_avg = econs_avg + md_ener%delta_cons!(md_ener%constant-constant)/nat*kelvin
       END SELECT

       ! Print to the screen
       iw = cp_print_key_unit_nr(logger,motion_section,"MD%PRINT%PROGRAM_RUN_INFO",&
            extension=".mdLog",file_position=my_pos, file_action=my_act,error=error)
       IF ( iw>0  ) THEN
          WRITE ( iw, '( )' )
          WRITE ( iw, '( A,A )' ) ' **************************************', &
               '*****************************************'
          WRITE ( iw, '( A,T61,A )' ) ' ENSEMBLE TYPE                = ', &
               ADJUSTR ( ensemble_name_types(simpar%ensemble)(1:20))
          WRITE ( iw, '( A,T71,I10 )' ) ' STEP NUMBER                  = ', itimes
          WRITE ( iw, '( A,T60,1(1X,F20.6) )' ) ' TIME [FS]                    = ', time*femtoseconds
          WRITE ( iw, '( A,T60,1(1X,E20.12) )' ) ' CONSERVED QNTY               = ', &
               md_ener%constant
          WRITE ( iw, '( )' )
          WRITE ( iw, '( T42,A )' ) 'INSTANTANEOUS        AVERAGES'
          WRITE ( iw, '( A,T39,2(1X,F20.2) )' ) ' CPU [S]                      = ', &
               used_time, averages%avecpu
          WRITE ( iw, '( A,T39,2(1X,E20.12) )' ) ' {E-E0}/{k_b*N_at}            = ', &
               econs, econs_avg/REAL( delta_t,dp)
          WRITE ( iw, '( A,A,T31,A,T39,2(1X,E20.12) )' ) ' POTENTIAL ENERGY', &
               ADJUSTL(unit % e_label), '= ', md_ener%epot, averages%avepot
          IF(simpar%ensemble/=reftraj_ensemble) THEN
            IF (.NOT.qmmm) THEN
              ! No QM/MM info
              WRITE ( iw, '( A,A,T31,A,T39,2(1X,E20.12) )' ) ' KINETIC ENERGY', &
                  ADJUSTL(unit % e_label), '= ', md_ener%ekin, averages%avekin
              WRITE ( iw, '( A,A,T31,A,T39,2(1X,F20.3) )' ) ' TEMPERATURE', &
                  ADJUSTL(unit % temp_label), '= ', md_ener%temp_part, averages%avetemp
            ELSE
              WRITE ( iw, '( A,A,T31,A,T39,2(1X,E20.12) )' ) ' MM KINETIC ENERGY', &
                  ADJUSTL(unit % e_label), '= ', md_ener%ekin, averages%avekin
              WRITE ( iw, '( A,A,T31,A,T39,2(1X,E20.12) )' ) ' QM KINETIC ENERGY', &
                  ADJUSTL(unit % e_label), '= ', md_ener%ekin_qm, averages%avekin_qm
              WRITE ( iw, '( A,A,T31,A,T39,2(1X,F20.3) )' ) ' MM TEMPERATURE', &
                  ADJUSTL(unit % temp_label), '= ', md_ener%temp_part, averages%avetemp
              WRITE ( iw, '( A,A,T31,A,T39,2(1X,F20.3) )' ) ' QM TEMPERATURE', &
                  ADJUSTL(unit % temp_label), '= ', md_ener%temp_qm, averages%avetemp_qm
            END IF
          END IF
          IF ( virial % pv_availability ) &
               WRITE ( iw, '( A,A,T31,A,T39,2(1X,E20.12) )' ) ' PRESSURE', &
               ADJUSTL(unit % pres_label), '= ', pv_scalar, averages%avepress
          IF ( simpar % ensemble == nph_ensemble .OR.&
               simpar % ensemble == nph_uniaxial_ensemble .OR.&
               simpar % ensemble == nph_uniaxial_damped_ensemble) THEN
             WRITE ( iw, '( A,A,T31,A,T39,2(1X,E20.12) )' ) ' P_XX', &
                  ADJUSTL(unit % pres_label), '= ', pv_xx, averages%avepxx
             WRITE ( iw, '( A,A,T31,A,T39,2(1X,E20.12) )' ) ' HUGONIOT', &
                  ADJUSTL(unit % temp_label), '= ', hugoniot/3._dp/nat*kelvin,&
                  averages%avehugoniot/3._dp/nat*kelvin
          END IF
          IF ( simpar % ensemble == nph_ensemble .OR.&
               simpar % ensemble == nph_uniaxial_ensemble .OR.&
               simpar % ensemble == nph_uniaxial_damped_ensemble.OR.&
               simpar % ensemble == npt_i_ensemble.OR.&
               simpar % ensemble == npt_f_ensemble ) THEN
             WRITE ( iw, '( A,A,T31,A,T39,2(1X,E20.12) )' ) ' BAROSTAT TEMP', &
                  ADJUSTL(unit % temp_label), '= ', md_ener%temp_baro, averages%avetemp_baro
             WRITE ( iw, '( A,A,T31,A,T39,2(1X,E20.12) )' ) ' VOLUME', &
                  ADJUSTL(unit % vol_label), '= ', cell%deth, averages%avevol
             WRITE ( iw, '( A,A,T31,A,T33,3(1X,E15.7) )' ) ' CELL LNTHS', &
                  ADJUSTL(unit % l_label), '= ', abc(1), abc(2), abc(3)
             WRITE ( iw, '( A,A,T31,A,T33,3(1X,E15.7) )' ) ' AVE. CELL LNTHS', &
                  ADJUSTL(unit % l_label), '= ', averages%aveca, averages%avecb, &
                  averages%avecc
          END IF
          IF (simpar%ensemble==npt_f_ensemble) THEN
             WRITE ( iw, '( A,A,T31,A,T33,3(1X,E15.7) )' ) ' CELL ANGLS', &
                  ADJUSTL(unit % angl_label), '= ', cell_angle(3), cell_angle(2), cell_angle(1)
             WRITE ( iw, '( A,A,T31,A,T33,3(1X,E15.7) )' ) ' AVE. CELL ANGLS', &
                  ADJUSTL(unit % angl_label), '= ', averages%aveal, averages%avebe, &
                  averages%avega
          END IF
          WRITE ( iw, '( A,A )' ) ' **************************************', &
               '*****************************************'
          WRITE ( iw, '( )' )
       END IF
       CALL cp_print_key_finished_output(iw,logger,motion_section,&
            "MD%PRINT%PROGRAM_RUN_INFO", error=error)
    END IF

    IF ((.NOT.init).OR.(itimes==0)) THEN
       ene = cp_print_key_unit_nr(logger,motion_section,"MD%PRINT%ENERGY",&
            extension=".ener",file_position=my_pos, file_action=my_act,error=error)

       CALL get_output_format(motion_section,"PRINT%TRAJECTORY",my_form,my_ext,error)
       crd = cp_print_key_unit_nr(logger,motion_section,"PRINT%TRAJECTORY",&
            extension=my_ext, file_position=my_pos, file_action=my_act,&
            file_form=my_form, middle_name="pos", is_new_file=init_traj,&
            error=error)
       IF (is_mixed) THEN
           ener_mix=cp_print_key_unit_nr(logger,motion_section,"PRINT%MIXED_ENERGIES",&
                extension=".ener", file_position=my_pos, file_action=my_act,&
                middle_name="mix", error=error)
       ELSE
           ener_mix=-1
       ENDIF
       IF    ( simpar % ensemble == nph_ensemble .OR.&
               simpar % ensemble == nph_uniaxial_ensemble .OR.&
               simpar % ensemble == nph_uniaxial_damped_ensemble.OR.&
               simpar % ensemble == npt_i_ensemble.OR.&
               simpar % ensemble == npt_f_ensemble ) THEN
         box = cp_print_key_unit_nr(logger,motion_section,"PRINT%CELL",&
              extension=".cell", file_position=my_pos, file_action=my_act,&
              file_form="FORMATTED", error=error)

       ELSE
         box = 0
       END IF
       CALL get_output_format(motion_section,"PRINT%VELOCITIES",my_form,my_ext,error)
       vel = cp_print_key_unit_nr(logger,motion_section,"PRINT%VELOCITIES",&
            extension=my_ext, file_position=my_pos, file_action=my_act,&
            file_form=my_form, middle_name="vel", is_new_file=init_vel,&
            error=error)
       CALL get_output_format(motion_section,"PRINT%FORCES",my_form,my_ext,error)
       frc = cp_print_key_unit_nr(logger,motion_section,"PRINT%FORCES",&
            extension=my_ext, file_position=my_pos, file_action=my_act,&
            file_form=my_form, middle_name="force", is_new_file=init_force,&
            error=error)
       IF( virial % pv_availability ) THEN
         str = cp_print_key_unit_nr(logger,motion_section,"PRINT%STRESS",&
              extension=".stress", file_position=my_pos, file_action=my_act,&
              file_form="FORMATTED", error=error)
       ELSE
         str = 0
       END IF
       IF(shell_present) THEN
         CALL force_env_get(force_env,harmonic_shell=harm_shell,error=error)
         shene=cp_print_key_unit_nr(logger,motion_section,"MD%PRINT%SHELL_ENERGY",&
            extension=".shener",file_position=my_pos, file_action=my_act,error=error)

         CALL get_output_format(motion_section,"PRINT%SHELL_TRAJECTORY",my_form,my_ext,error)
         shcrd =  cp_print_key_unit_nr(logger,motion_section,"PRINT%SHELL_TRAJECTORY",&
            extension=my_ext, file_position=my_pos, file_action=my_act,&
            file_form=my_form, middle_name="shpos", is_new_file=init_shtraj,&
            error=error)

         CALL cp_subsys_get(subsys,shell_particles=shell_particles,error=error)
         nshell=shell_particles%n_els
       ELSE
         shene = 0
         shcrd = 0
       END IF
       IF(shell_adiabatic) THEN
         CALL get_output_format(motion_section,"PRINT%CORE_TRAJECTORY",my_form,my_ext,error)
         cocrd = cp_print_key_unit_nr(logger,motion_section,"PRINT%CORE_TRAJECTORY",&
            extension=my_ext, file_position=my_pos, file_action=my_act,&
            file_form=my_form, middle_name="copos", is_new_file=init_cotraj,&
            error=error)

         CALL get_output_format(motion_section,"PRINT%SHELL_VELOCITIES",my_form,my_ext,error)
         shvel = cp_print_key_unit_nr(logger,motion_section,"PRINT%SHELL_VELOCITIES",&
            extension=my_ext, file_position=my_pos, file_action=my_act,&
            file_form=my_form, middle_name="shvel", is_new_file=init_shvel,&
            error=error)
         
         CALL get_output_format(motion_section,"PRINT%CORE_VELOCITIES",my_form,my_ext,error)
         covel = cp_print_key_unit_nr(logger,motion_section,"PRINT%CORE_VELOCITIES",&
            extension=my_ext, file_position=my_pos, file_action=my_act,&
            file_form=my_form, middle_name="covel", is_new_file=init_covel,&
            error=error)

         CALL get_output_format(motion_section,"PRINT%SHELL_FORCES",my_form,my_ext,error)
         shfrc = cp_print_key_unit_nr(logger,motion_section,"PRINT%SHELL_FORCES",&
            extension=my_ext, file_position=my_pos, file_action=my_act,&
            file_form=my_form, middle_name="shfrc", is_new_file=init_shforce,&
            error=error)

         CALL get_output_format(motion_section,"PRINT%CORE_FORCES",my_form,my_ext,error)
         cofrc = cp_print_key_unit_nr(logger,motion_section,"PRINT%CORE_FORCES",&
            extension=my_ext, file_position=my_pos, file_action=my_act,&
            file_form=my_form, middle_name="cofrc", is_new_file=init_coforce,&
            error=error)

       ELSE
         cocrd = 0 
         shvel = 0 
         covel = 0
         shfrc = 0 
         cofrc = 0
       END IF
       IF(simpar%temperature_per_kind) THEN
          tempkind=cp_print_key_unit_nr(logger,motion_section,"MD%PRINT%TEMP_KIND",&
                   extension=".temp",file_position=my_pos, file_action=my_act,error=error)
          tempshellkind = 0
          IF(shell_adiabatic) THEN
             tempshellkind=cp_print_key_unit_nr(logger,motion_section,&
                           "MD%PRINT%TEMP_SHELL_KIND",extension=".shtemp",&
                           file_position=my_pos, file_action=my_act,error=error)
          END IF
       ELSE
          tempkind = 0
          tempshellkind = 0
       END IF
       
       !-------------------------------------------
       ! Print QMMM translation vector if requested
       !-------------------------------------------
       IF (qmmm) THEN
          trsl = cp_print_key_unit_nr(logger,motion_section,"PRINT%TRANSLATION_VECTOR",&
               extension=".translation", middle_name="qmmm", error=error)
          IF (trsl>0) THEN
             WRITE(trsl,'(I10,3F15.10)')itimes,force_env%qmmm_env%transl_v
          END IF
          CALL cp_print_key_finished_output(trsl,logger,motion_section,&
               "PRINT%TRANSLATION_VECTOR", error=error)
       END IF
       !------------------
       ! Print Coordinates
       !------------------
       CALL write_structure_data(particles%els,cell,motion_section,error)
       IF ( crd > 0 ) THEN
         CALL section_vals_val_get(motion_section,"PRINT%TRAJECTORY%FORMAT",i_val=outformat,error=error)
          IF(outformat==dump_dcd) THEN
             IF ( init_traj ) THEN
                !Lets write the header for the coordinate dcd
                iskip = section_get_ival(motion_section,"PRINT%TRAJECTORY%EACH",error=error)
                WRITE(crd) "CORD",0,itimes,iskip,&
                     0,0,0,0,0,0,&
                     REAL(simpar%dt*femtoseconds,KIND=sp),1,0,0,0,0,0,0,0,0,24
                remark1= "REMARK FILETYPE CORD DCD GENERATED BY CP2K"
                remark2= "REMARK Support new DCD format with cell information"
                WRITE(crd) 2,remark1,remark2
                WRITE(crd) nat
                CALL m_flush(crd)
             ENDIF
          ELSEIF(outformat==dump_xmol) THEN
             !         *** This information can be digested by Molden ***
             WRITE (UNIT=title,FMT="(A,I8,A,F20.10,A,F12.3)")&
                  " i = ",itimes," , E =",md_ener%epot,", Time (fs)=",time*femtoseconds
          ELSE
             title = ""
          END IF
          CALL cp_subsys_get(subsys,particles=particles,error=error)
          CALL write_particle_coordinates(particles%els,crd,&
               outformat,"POS",title,cell)
          CALL m_flush(crd)
       END IF
       !------------------
       ! Print Simulation Cell (only if constant pressure)
       !------------------
       IF ( box > 0 ) THEN
          WRITE(UNIT=box,FMT="(I8, F12.3,10F20.10)") itimes,time*femtoseconds,&
             cell%hmat(1,1)*angstrom,cell%hmat(2,1)*angstrom,cell%hmat(3,1)*angstrom,&
             cell%hmat(1,2)*angstrom,cell%hmat(2,2)*angstrom,cell%hmat(3,2)*angstrom,&
             cell%hmat(1,3)*angstrom,cell%hmat(2,3)*angstrom,cell%hmat(3,3)*angstrom,&
             cell%deth*angstrom*angstrom*angstrom
      
          CALL m_flush(box)
       END IF
       !------------------
       ! Print the Stress tensor
       !------------------
       IF ( str > 0 ) THEN
          WRITE(UNIT=str,FMT="(I10,F12.3,9F20.10)") itimes,time*femtoseconds,&
                        virial%pv_total(1,1)/cell%deth*unit%bar,&
                        virial%pv_total(1,2)/cell%deth*unit%bar,&
                        virial%pv_total(1,3)/cell%deth*unit%bar,&
                        virial%pv_total(2,1)/cell%deth*unit%bar,&
                        virial%pv_total(2,2)/cell%deth*unit%bar,&
                        virial%pv_total(2,3)/cell%deth*unit%bar,& 
                        virial%pv_total(3,1)/cell%deth*unit%bar,&
                        virial%pv_total(3,2)/cell%deth*unit%bar,&
                        virial%pv_total(3,3)/cell%deth*unit%bar
          CALL m_flush(str)
       END IF
       !------------------
       ! Print Velocities
       !------------------
       IF ( vel > 0 ) THEN
          CALL section_vals_val_get(motion_section,"PRINT%VELOCITIES%FORMAT",i_val=outformat,error=error)
          IF(outformat==dump_dcd) THEN
             IF ( init_vel ) THEN
                iskip = section_get_ival(motion_section,"PRINT%VELOCITIES%EACH",error=error)
                WRITE(vel) "VEL ",0,itimes,iskip,&
                     0,0,0,0,0,0,&
                     REAL(simpar%dt*femtoseconds,KIND=sp), 1,0,0,0,0,0,0,0,0,24
                remark1= "REMARK FILETYPE VEL  DCD GENERATED BY CP2K"
                remark2= "REMARK Support new DCD format with cell information"
                WRITE(vel) 2,remark1,remark2
                WRITE(vel) nat
                CALL m_flush(vel)
             END IF
          ELSEIF(outformat==dump_xmol) THEN
             !         *** This information can be digested by Molden ***
             WRITE (UNIT=title,FMT="(A,I8,A,F20.10)")&
                  " i =",itimes,", E =",md_ener%epot
          ELSE
             title = ""
          END IF
          CALL cp_subsys_get(subsys,particles=particles,error=error)
          CALL write_particle_coordinates(particles%els,vel,&
               outformat,"VEL",title,cell)
          CALL m_flush(vel)
       END IF
       !------------------
       ! Print Force
       !------------------
       IF ( frc > 0 ) THEN
          CALL section_vals_val_get(motion_section,"PRINT%FORCES%FORMAT",i_val=outformat,error=error)
          IF(outformat==dump_dcd) THEN
             IF ( init_force ) THEN
                iskip = section_get_ival(motion_section,"PRINT%FORCES%EACH",error=error)
                WRITE(frc) "FRC ",0,itimes,iskip,&
                     0,0,0,0,0,0,&
                     REAL(simpar%dt*femtoseconds,KIND=sp), 1,0,0,0,0,0,0,0,0,24
                remark1= "REMARK FILETYPE FRC  DCD GENERATED BY CP2K"
                remark2= "REMARK Support new DCD format with cell information"
                WRITE(frc) 2,remark1,remark2
                WRITE(frc) nat
                CALL m_flush(frc)
             END IF
          ELSEIF(outformat==dump_xmol) THEN
             !         *** This information can be digested by Molden ***
             WRITE (UNIT=title,FMT="(A,I8,A,F20.10)")&
                  " i =",itimes,", E =",md_ener%epot
          ELSE
             title = ""
          END IF
          CALL cp_subsys_get(subsys,particles=particles,error=error)
          CALL write_particle_coordinates(particles%els,frc,&
               outformat,"FORCE",title,cell)
          CALL m_flush(frc)
       END IF
       !------------------
       ! Print MIXED Energy
       !------------------
       IF (ener_mix>0 .AND. is_mixed) THEN
          WRITE (ener_mix,"(I8,F12.3,F20.9,"//&
               cp_to_string(SIZE(force_env%mixed_env%energies))//"F20.9,F20.9)")&
               itimes,time*femtoseconds,md_ener%epot,force_env%mixed_env%energies,md_ener%constant
          CALL m_flush(ener_mix)
       END IF
       !------------------
       ! Print Energy
       !------------------
       IF (ene>0) THEN
          ! please change also the corresponding format explaination above
          ! keep the constant of motion the true constant of motion !
          ! "Step Nr.    Time[fs]  Kin.[a.u.]     Temp[K]  Pot.[a.u.]    Cons Qty UsedTime[s]"
          !
          WRITE (ene,"(I10,F20.3,F20.9,F20.9,F20.9,F20.9,F20.9)")&
               itimes,time*femtoseconds,md_ener%ekin,md_ener%temp_part, md_ener%epot,md_ener%constant,used_time
          CALL m_flush(ene)
       END IF
       !------------------
       ! Print Shell Energy
       !------------------
       IF (shene>0) THEN
          ! "Step Nr.    Time[fs]  Shell Kin.[a.u.]    Shell Temp[K]  Shell Pot.[a.u.] 
          !
          WRITE (shene,"(I8,F12.3,F20.9,F20.9,F20.9,F20.9 )")&
               itimes,time*femtoseconds,md_ener%ekin_shell,md_ener%temp_shell,harm_shell
          CALL m_flush(shene)
       END IF
       !------------------
       !  Some additional printing for the shell-model (can be displaced in another routine)
       !------------------
       ! Print Shell Coordinates
       !------------------
       IF ( shcrd > 0 ) THEN
          CALL section_vals_val_get(motion_section,"PRINT%SHELL_TRAJECTORY%FORMAT",i_val=outformat,error=error)
          IF(outformat==dump_dcd) THEN
             IF ( init_shtraj ) THEN
                !Lets write the header for the coordinate dcd
                iskip = section_get_ival(motion_section,"PRINT%SHELL_TRAJECTORY%EACH",error=error)
                WRITE(shcrd) "CORD",0,itimes,iskip,&
                     0,0,0,0,0,0,&
                     REAL(simpar%dt*femtoseconds,KIND=sp),1,0,0,0,0,0,0,0,0,24
                remark1= "REMARK FILETYPE CORD DCD GENERATED BY CP2K"
                remark2= "REMARK Support new DCD format with cell information"
                WRITE(shcrd) 2,remark1,remark2
                WRITE(shcrd) nshell
                CALL m_flush(shcrd)
             ENDIF
          ELSEIF(outformat==dump_xmol) THEN
             !         *** This information can be digested by Molden ***
             WRITE (UNIT=title,FMT="(A,I8,A,F20.10)")&
                  " i =",itimes,", E =",md_ener%epot
          ELSE
             title = ""
          END IF
          CALL cp_subsys_get(subsys,shell_particles=shell_particles,error=error)
          CALL write_particle_coordinates(shell_particles%els,shcrd,&
               outformat,"POS",title,cell)
          CALL m_flush(shcrd)
       END IF
       !------------------
       ! Print Shell Velocities
       !------------------
       IF ( shvel > 0 ) THEN
          CALL section_vals_val_get(motion_section,"PRINT%SHELL_VELOCITIES%FORMAT",i_val=outformat,error=error)
          IF(outformat==dump_dcd) THEN
             IF ( init_shvel ) THEN
                iskip = section_get_ival(motion_section,"PRINT%SHELL_VELOCITIES%EACH",error=error)
                WRITE(vel) "VEL ",0,itimes,iskip,&
                     0,0,0,0,0,0,&
                     REAL(simpar%dt*femtoseconds,KIND=sp), 1,0,0,0,0,0,0,0,0,24
                remark1= "REMARK FILETYPE VEL  DCD GENERATED BY CP2K"
                remark2= "REMARK Support new DCD format with cell information"
                WRITE(shvel) 2,remark1,remark2
                WRITE(shvel) nshell
                CALL m_flush(shvel)
             END IF
          ELSEIF(outformat==dump_xmol) THEN
             !         *** This information can be digested by Molden ***
             WRITE (UNIT=title,FMT="(A,I8,A,F20.10)")&
                  " i =",itimes,", E =",md_ener%epot
          ELSE
             title = ""
          END IF
          CALL cp_subsys_get(subsys,shell_particles=shell_particles,error=error)
          CALL write_particle_coordinates(shell_particles%els,shvel,&
               outformat,"VEL",title,cell)
          CALL m_flush(shvel)
       END IF
       !------------------
       ! Print Shell Force
       !------------------
       IF ( shfrc > 0 ) THEN
          CALL section_vals_val_get(motion_section,"PRINT%SHELL_FORCES%FORMAT",i_val=outformat,error=error)
          IF(outformat==dump_dcd) THEN
             IF ( init_shforce ) THEN
                iskip = section_get_ival(motion_section,"PRINT%SHELL_FORCES%EACH",error=error)
                WRITE(shfrc) "FRC ",0,itimes,iskip,&
                     0,0,0,0,0,0,&
                     REAL(simpar%dt*femtoseconds,KIND=sp), 1,0,0,0,0,0,0,0,0,24
                remark1= "REMARK FILETYPE FRC  DCD GENERATED BY CP2K"
                remark2= "REMARK Support new DCD format with cell information"
                WRITE(shfrc) 2,remark1,remark2
                WRITE(shfrc) nshell
                CALL m_flush(shfrc)
             END IF
          ELSEIF(outformat==dump_xmol) THEN
             !         *** This information can be digested by Molden ***
             WRITE (UNIT=title,FMT="(A,I8,A,F20.10)")&
                  " i =",itimes,", E =",md_ener%epot
          ELSE
             title = ""
          END IF
          CALL cp_subsys_get(subsys,shell_particles=shell_particles,error=error)
          CALL write_particle_coordinates(shell_particles%els,shfrc,&
               outformat,"FORCE",title,cell)
          CALL m_flush(shfrc)
       END IF
       !------------------
       ! Print Core Coordinates
       !------------------
       IF ( cocrd > 0 ) THEN
          CALL section_vals_val_get(motion_section,"PRINT%CORE_TRAJECTORY%FORMAT",i_val=outformat,error=error)
          IF(outformat==dump_dcd) THEN
             IF ( init_cotraj ) THEN
                !Lets write the header for the coordinate dcd
                iskip = section_get_ival(motion_section,"PRINT%CORE_TRAJECTORY%EACH",error=error)
                WRITE(cocrd) "CORD",0,itimes,iskip,&
                     0,0,0,0,0,0,&
                     REAL(simpar%dt*femtoseconds,KIND=sp),1,0,0,0,0,0,0,0,0,24
                remark1= "REMARK FILETYPE CORD DCD GENERATED BY CP2K"
                remark2= "REMARK Support new DCD format with cell information"
                WRITE(cocrd) 2,remark1,remark2
                WRITE(cocrd) nshell
                CALL m_flush(cocrd)
             ENDIF
          ELSEIF(outformat==dump_xmol) THEN
             !         *** This information can be digested by Molden ***
             WRITE (UNIT=title,FMT="(A,I8,A,F20.10)")&
                  " i =",itimes,", E =",md_ener%epot
          ELSE
             title = ""
          END IF
          CALL cp_subsys_get(subsys,core_particles=core_particles,error=error)
          CALL write_particle_coordinates(core_particles%els,cocrd,&
               outformat,"POS",title,cell)
          CALL m_flush(cocrd)
       END IF
       !------------------
       ! Print core Velocities
       !------------------
       IF ( covel > 0 ) THEN
          CALL section_vals_val_get(motion_section,"PRINT%CORE_VELOCITIES%FORMAT",i_val=outformat,error=error)
          IF(outformat==dump_dcd) THEN
             IF ( init_covel ) THEN
                iskip = section_get_ival(motion_section,"PRINT%CORE_VELOCITIES%EACH",error=error)
                WRITE(vel) "VEL ",0,itimes,iskip,&
                     0,0,0,0,0,0,&
                     REAL(simpar%dt*femtoseconds,KIND=sp), 1,0,0,0,0,0,0,0,0,24
                remark1= "REMARK FILETYPE VEL  DCD GENERATED BY CP2K"
                remark2= "REMARK Support new DCD format with cell information"
                WRITE(covel) 2,remark1,remark2
                WRITE(covel) nshell
                CALL m_flush(covel)
             END IF
          ELSEIF(outformat==dump_xmol) THEN
             !         *** This information can be digested by Molden ***
             WRITE (UNIT=title,FMT="(A,I8,A,F20.10)")&
                  " i =",itimes,", E =",md_ener%epot
          ELSE
             title = ""
          END IF
          CALL cp_subsys_get(subsys,core_particles=core_particles,error=error)
          CALL write_particle_coordinates(core_particles%els,covel,&
               outformat,"VEL",title,cell)
          CALL m_flush(covel)
       END IF
       !------------------
       ! Print core Force
       !------------------
       IF ( cofrc > 0 ) THEN
          CALL section_vals_val_get(motion_section,"PRINT%CORE_FORCES%FORMAT",i_val=outformat,error=error)
          IF(outformat==dump_dcd) THEN
             IF ( init_coforce ) THEN
                iskip = section_get_ival(motion_section,"PRINT%CORE_FORCES%EACH",error=error)
                WRITE(cofrc) "FRC ",0,itimes,iskip,&
                     0,0,0,0,0,0,&
                     REAL(simpar%dt*femtoseconds,KIND=sp), 1,0,0,0,0,0,0,0,0,24
                remark1= "REMARK FILETYPE FRC  DCD GENERATED BY CP2K"
                remark2= "REMARK Support new DCD format with cell information"
                WRITE(cofrc) 2,remark1,remark2
                WRITE(cofrc) nshell
                CALL m_flush(cofrc)
             END IF
          ELSEIF(outformat==dump_xmol) THEN
             !         *** This information can be digested by Molden ***
             WRITE (UNIT=title,FMT="(A,I8,A,F20.10)")&
                  " i =",itimes,", E =",md_ener%epot
          ELSE
             title = ""
          END IF
          CALL cp_subsys_get(subsys,core_particles=core_particles,error=error)
          CALL write_particle_coordinates(core_particles%els,cofrc,&
               outformat,"FORCE",title,cell)
          CALL m_flush(cofrc)
       END IF

       IF( tempkind > 0 ) THEN
          nkind = SIZE(md_ener%temp_kind)
          fmd="(I10,F20.3,"//TRIM(ADJUSTL(cp_to_string(nkind)))//"F20.9)"     
          fmd=TRIM(fmd)
          WRITE (tempkind,fmd)&
               itimes,time*femtoseconds, md_ener%temp_kind(1:nkind)
          CALL m_flush(tempkind)
       END IF   
       IF( tempshellkind > 0 ) THEN
          nkind = SIZE(md_ener%temp_shell_kind)
          fmd="(I10,F20.3,"//TRIM(ADJUSTL(cp_to_string(nkind)))//"F20.9)"     
          fmd=TRIM(fmd)
          WRITE (tempshellkind,fmd)&
               itimes,time*femtoseconds, md_ener%temp_shell_kind(1:nkind)
          CALL m_flush(tempshellkind)
       END IF

       CALL cp_print_key_finished_output(ene,logger,motion_section,&
            "MD%PRINT%ENERGY", error=error)
       CALL cp_print_key_finished_output(crd,logger,motion_section,&
            "PRINT%TRAJECTORY", error=error)
       IF ( simpar % ensemble == nph_ensemble .OR.&
               simpar % ensemble == nph_uniaxial_ensemble .OR.&
               simpar % ensemble == nph_uniaxial_damped_ensemble.OR.&
               simpar % ensemble == npt_i_ensemble.OR.&
               simpar % ensemble == npt_f_ensemble ) THEN
         CALL cp_print_key_finished_output(box,logger,motion_section,&
              "PRINT%CELL", error=error)
       END IF
       CALL cp_print_key_finished_output(vel,logger,motion_section,&
            "PRINT%VELOCITIES", error=error)
       CALL cp_print_key_finished_output(frc,logger,motion_section,&
            "PRINT%FORCES", error=error)
       IF (is_mixed) THEN
           CALL cp_print_key_finished_output(ener_mix,logger,motion_section,&
                "PRINT%MIXED_ENERGIES", error=error)
       ENDIF
       IF( virial % pv_availability ) THEN
         CALL cp_print_key_finished_output(str,logger,motion_section,&
              "PRINT%STRESS", error=error)
       END IF
       IF(shell_present) THEN
         CALL cp_print_key_finished_output(shene,logger,motion_section,&
            "MD%PRINT%SHELL_ENERGY", error=error)
         CALL cp_print_key_finished_output(shcrd,logger,motion_section,&
            "PRINT%SHELL_TRAJECTORY", error=error)
       END IF
       IF(shell_adiabatic) THEN
         CALL cp_print_key_finished_output(cocrd,logger,motion_section,&
            "PRINT%CORE_TRAJECTORY", error=error)
         CALL cp_print_key_finished_output(shvel,logger,motion_section,&
            "PRINT%SHELL_VELOCITIES", error=error)
         CALL cp_print_key_finished_output(covel,logger,motion_section,&
            "PRINT%CORE_VELOCITIES", error=error)
         CALL cp_print_key_finished_output(shfrc,logger,motion_section,&
            "PRINT%SHELL_FORCES", error=error)
         CALL cp_print_key_finished_output(cofrc,logger,motion_section,&
            "PRINT%CORE_FORCES", error=error)
       END IF
       IF(simpar%temperature_per_kind) THEN
         CALL cp_print_key_finished_output(tempkind,logger,motion_section,&
              "MD%PRINT%TEMP_KIND", error=error) 
         IF(shell_adiabatic) THEN
           CALL cp_print_key_finished_output(tempshellkind,logger,motion_section,&
              "MD%PRINT%TEMP_SHELL_KIND", error=error) 
         END IF
       END IF
    END IF  ! not init
    init = .FALSE.
    CALL set_md_env(md_env,init=init,error=error)

    CALL timestop(handle)

  END SUBROUTINE md_write_output

!******************************************************************************
!!****** md_energies/get_output_format [1.0] *
!!
!!   NAME
!!     get_output_format
!!
!!   FUNCTION
!!     Info on the unit to be opened to dump MD informations
!!
!!   AUTHOR
!!     Teodoro Laino - University of Zurich - 07.2007
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*****************************************************************************
  SUBROUTINE get_output_format(section,path,my_form,my_ext,error)
    TYPE(section_vals_type), POINTER         :: section
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: path
    CHARACTER(LEN=*), INTENT(OUT)            :: my_form, my_ext
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_output_format', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: output_format

    my_form = "FORMATTED"
    my_ext  = ".xyz"
    IF (PRESENT(path)) THEN
       CALL section_vals_val_get(section,TRIM(path)//"%FORMAT",i_val=output_format,error=error)
    ELSE
       CALL section_vals_val_get(section,"FORMAT",i_val=output_format,error=error)
    END IF
    IF(output_format==dump_dcd) THEN
       my_form = "UNFORMATTED"
       my_ext  = ".dcd"
    END IF
  END SUBROUTINE get_output_format
 
END MODULE md_energies
