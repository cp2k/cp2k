!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief prints all energy info per timestep to the screen or to
!>      user defined output files
!> \author Joost VandeVondele (copy from md_fist_energies)
!> 
!> \par History
!>      - New MD data are appended to the old data (15.09.2003,MK)
! *****************************************************************************
MODULE md_energies
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE averages_types,                  ONLY: average_quantities_type,&
                                             compute_averages
  USE barostat_types,                  ONLY: barostat_type
  USE barostat_utils,                  ONLY: print_barostat_status
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell
  USE cp2k_units,                      ONLY: get_cp2k_units,&
                                             unit_convert_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_type
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type,&
                                             use_mixed_force
  USE input_constants,                 ONLY: ensemble_name_types,&
                                             npe_f_ensemble,&
                                             npe_i_ensemble,&
                                             nph_uniaxial_damped_ensemble,&
                                             nph_uniaxial_ensemble,&
                                             npt_f_ensemble,&
                                             npt_i_ensemble,&
                                             reftraj_ensemble
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE machine,                         ONLY: m_flush
  USE md_conserved_quantities,         ONLY: compute_conserved_quantity
  USE md_ener_types,                   ONLY: md_ener_type,&
                                             zero_md_ener
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type,&
                                             set_md_env
  USE motion_utils,                    ONLY: write_simulation_cell,&
                                             write_stress_tensor,&
                                             write_trajectory
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: write_structure_data
  USE physcon,                         ONLY: angstrom,&
                                             femtoseconds,&
                                             kelvin
  USE qmmm_types,                      ONLY: qmmm_env_qm_type
  USE reftraj_types,                   ONLY: reftraj_type
  USE simpar_types,                    ONLY: simpar_type
  USE thermostat_types,                ONLY: thermostats_type
  USE thermostat_utils,                ONLY: print_thermostats_status
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'md_energies'


  PUBLIC :: initialize_md_ener,&
            md_energy,&
            md_ener_reftraj,&
            md_write_output

CONTAINS

! *****************************************************************************
!> \par History
!>      -10-2007 created 
!> \author MI
! *****************************************************************************
  SUBROUTINE initialize_md_ener(md_ener,force_env,simpar,error)

    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(simpar_type), POINTER               :: simpar
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

! *****************************************************************************
!> \par History
!>      -10-2007 created 
!> \author MI
! *****************************************************************************
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
    TYPE(simpar_type), POINTER               :: simpar

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
       IF ( &
            simpar%ensemble == nph_uniaxial_ensemble .OR.&
            simpar%ensemble == nph_uniaxial_damped_ensemble ) THEN
          IF ( simpar%v0 == 0.0_dp ) simpar%v0 = cell_ref%deth
       ENDIF
    ENDIF

    CALL compute_conserved_quantity(md_env, md_ener, tkind=tkind,&
         tcoeff=tcoeff, tshell=tshell, natom=natom, error=error)
    
  END SUBROUTINE md_energy

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created 
!> \author MI
! *****************************************************************************
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
       md_ener%delta_epot = (reftraj%epot - reftraj%epot0)/REAL(reftraj%natom, kind=dp)*kelvin
    END IF

  END SUBROUTINE  md_ener_reftraj

! *****************************************************************************
!> \brief This routine computes the conserved quantity, temperature
!>      and things like that and prints them out
!> \par History
!>      - New MD data are appended to the old data (15.09.2003,MK)
!>      - 02.2008 - Teodoro Laino [tlaino] - University of Zurich
!>                  Cleaning code and collecting the many commons routines..
!> \author CJM
! *****************************************************************************
  SUBROUTINE md_write_output(md_env, error)

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'md_write_output', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: fmd, my_act, my_pos
    INTEGER                                  :: delta_t, ene, ener_mix, &
                                                handle, i, nat, nkind, shene, &
                                                tempkind, trsl
    INTEGER, POINTER                         :: itimes
    LOGICAL                                  :: failure, init, is_mixed, &
                                                new_file, qmmm, &
                                                shell_adiabatic, shell_present
    REAL(dp)                                 :: abc( 3 ), cell_angle( 3 ), &
                                                dt, econs, pv_scalar, pv_xx
    REAL(KIND=dp)                            :: harm_shell, hugoniot
    REAL(KIND=dp), POINTER                   :: time, used_time
    REAL(KIND=dp), SAVE                      :: itimes_old, t_old
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(average_quantities_type), POINTER   :: averages
    TYPE(barostat_type), POINTER             :: barostat
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(reftraj_type), POINTER              :: reftraj
    TYPE(section_vals_type), POINTER         :: motion_section, print_key, &
                                                root_section
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(thermostats_type), POINTER          :: thermostats
    TYPE(unit_convert_type), POINTER         :: unit
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.
    NULLIFY(logger,unit)
    logger => cp_error_get_logger(error)
    CALL timeset(routineN,handle)
    CALL get_cp2k_units ( unit )

    ! Zeroing
    hugoniot = 0.0_dp
    econs    = 0.0_dp
    shell_adiabatic = .FALSE.
    shell_present   = .FALSE.
    NULLIFY(motion_section, atomic_kinds, atomic_kind_set)
    NULLIFY(cell, subsys, force_env, md_ener, qmmm_env, reftraj)
    NULLIFY(core_particles, particles,shell_particles, print_key)
    NULLIFY(root_section,simpar,virial,thermostats)

    CALL get_md_env(md_env=md_env, itimes=itimes, t=time, used_time=used_time,&
         simpar=simpar, force_env=force_env, init=init, md_ener=md_ener,&
         reftraj=reftraj, thermostats=thermostats, barostat=barostat, &
         para_env=para_env, averages=averages, error=error)

    root_section   => force_env%root_section
    motion_section => section_vals_get_subs_vals(root_section,"MOTION",error=error)

    CALL force_env_get(force_env, cell=cell, subsys=subsys, qmmm_env=qmmm_env,&
         virial=virial, error=error)
    qmmm = .FALSE.
    IF (md_ener%nfree_qm>0) qmmm = .TRUE.

    is_mixed = (force_env%in_use == use_mixed_force)

    CALL cp_subsys_get(subsys,particles=particles,error=error)
    nat = particles%n_els
    dt  = simpar%dt*simpar%dt_fact

    ! Computing the scalar pressure
    IF ( virial%pv_availability ) THEN
       pv_scalar = 0._dp
       DO i = 1, 3
          pv_scalar = pv_scalar + virial%pv_total(i,i)
       END DO
       pv_scalar = pv_scalar/3._dp/cell%deth
       pv_scalar = pv_scalar * unit%bar
       pv_xx     = virial%pv_total(1,1)/cell%deth*unit%bar
    ENDIF

    CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds, error=error)
    atomic_kind_set => atomic_kinds%els
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                             shell_present=shell_present, &
                             shell_adiabatic=shell_adiabatic)

    CALL get_cell (cell, abc=abc, alpha=cell_angle(3), beta=cell_angle(2), gamma=cell_angle(1))

    ! Determine POS and ACT for I/O
    my_pos = "APPEND"
    my_act = "WRITE"
    IF (init.AND.(itimes==0)) THEN
       my_pos = "REWIND"
       my_act = "WRITE"
    END IF

    ! In case of REFTRAJ ensemble the POS is determined differently..
    ! according the REFTRAJ counter
    IF (simpar%ensemble==reftraj_ensemble) THEN
       IF((reftraj%isnap==reftraj%info%first_snapshot)) THEN
          my_pos  = "REWIND"
       END IF
    ENDIF

    ! Performing protocol relevant to the first step of an MD run
    IF (init) THEN
       ! Computing the Hugoniot for NPH calculations
       IF ( simpar%ensemble == nph_uniaxial_ensemble .OR.&
            simpar%ensemble == nph_uniaxial_damped_ensemble ) THEN
          IF (simpar%e0 == 0._dp ) simpar%e0 = md_ener%epot + md_ener%ekin
          hugoniot = md_ener%epot + md_ener%ekin - simpar%e0 - 0.5_dp*(pv_xx/unit%bar + simpar%p0)*&
               (simpar%v0-cell%deth)
       ENDIF

       ! Assigning the initial time-step and time
       itimes_old = itimes
       t_old      = time

       IF (simpar%ensemble==reftraj_ensemble) THEN
          itimes_old = reftraj%itimes0
          t_old = reftraj%time0
          reftraj%init = init
       END IF
    ELSE
       ! Performing protocol for anything beyond the first step of MD
       IF ( simpar%ensemble == nph_uniaxial_ensemble .OR. simpar%ensemble == nph_uniaxial_damped_ensemble) THEN
          hugoniot = md_ener%epot + md_ener%ekin - simpar%e0 - 0.5_dp * ( pv_xx/unit%bar + simpar%p0 ) * &
               ( simpar%v0 - cell%deth )
       END IF

       IF (simpar%ensemble==reftraj_ensemble) THEN
          time = reftraj%time
          econs = md_ener%delta_epot
          itimes = reftraj%itimes
       ELSE
          econs = md_ener%delta_cons
       END IF

       ! Compute average quantities
       delta_t = itimes - itimes_old
       CALL compute_averages(averages, force_env, md_ener, cell, virial, pv_scalar, pv_xx,&
            used_time, hugoniot, abc, cell_angle, delta_t, nat, error)
    END IF
    
    ! Print md information
    CALL md_write_info_low(simpar, md_ener, unit, qmmm, virial, reftraj, cell, abc,&
         cell_angle, itimes, dt, time, used_time, averages, econs, pv_scalar, pv_xx,&
         hugoniot, nat, init, logger, motion_section, my_pos, my_act, error)

    ! Real Ouput driven by the PRINT sections
    IF ((.NOT.init).OR.(itimes==0).OR.simpar%ensemble==reftraj_ensemble) THEN
       ! Print Energy
       ene = cp_print_key_unit_nr(logger,motion_section,"MD%PRINT%ENERGY",&
            extension=".ener",file_position=my_pos, file_action=my_act, is_new_file=new_file,&
            error=error)
       IF (ene>0) THEN
          IF (new_file) THEN
             ! Please change also the corresponding format explaination below 
             ! keep the constant of motion the true constant of motion !
             WRITE (ene,'("#",5X,A,10X,A,8X,A,10X,A,12X,A,2(8X,A))')"Step Nr.","Time[fs]","Kin.[a.u.]","Temp[K]",&
                  "Pot.[a.u.]","Cons Qty[a.u.]","UsedTime[s]"
          END IF
          WRITE (ene,"(I10,F20.3,F20.9,F20.9,F20.9,F20.9,F20.9)")&
               itimes,time*femtoseconds,md_ener%ekin,md_ener%temp_part, md_ener%epot,md_ener%constant,used_time
          CALL m_flush(ene)
       END IF
       CALL cp_print_key_finished_output(ene,logger,motion_section,"MD%PRINT%ENERGY", error=error)

       ! Possibly Print MIXED Energy 
       IF (is_mixed) THEN
          ener_mix=cp_print_key_unit_nr(logger,motion_section,"PRINT%MIXED_ENERGIES",&
               extension=".ener", file_position=my_pos, file_action=my_act,&
               middle_name="mix", error=error)
          IF (ener_mix>0) THEN
             WRITE (ener_mix,"(I8,F12.3,F20.9,"//cp_to_string(SIZE(force_env%mixed_env%energies))//"F20.9,F20.9)")&
                  itimes,time*femtoseconds,md_ener%epot,force_env%mixed_env%energies,md_ener%constant
             CALL m_flush(ener_mix)
          END IF
          CALL cp_print_key_finished_output(ener_mix,logger,motion_section,"PRINT%MIXED_ENERGIES", error=error)
       ENDIF

       ! Print QMMM translation vector if requested
       IF (qmmm) THEN
          trsl = cp_print_key_unit_nr(logger,motion_section,"PRINT%TRANSLATION_VECTOR",&
               extension=".translation", middle_name="qmmm", error=error)
          IF (trsl>0) THEN
             WRITE(trsl,'(I10,3F15.10)')itimes,qmmm_env%transl_v
          END IF
          CALL cp_print_key_finished_output(trsl,logger,motion_section,&
               "PRINT%TRANSLATION_VECTOR", error=error)
       END IF

       ! Write Structure data
       CALL write_structure_data(particles%els,cell,motion_section,error)

       ! Print Coordinates
       CALL write_trajectory(force_env, root_section, itimes, time*femtoseconds, dt*femtoseconds, md_ener%epot,&
            pos=my_pos, act=my_act, extended_xmol_title=.TRUE., error=error)

       ! Print Velocities
       CALL write_trajectory(force_env, root_section, itimes, time*femtoseconds, dt*femtoseconds, md_ener%epot,&
            "VELOCITIES", my_pos, my_act, middle_name="vel", extended_xmol_title=.TRUE., error=error)

       ! Print Force
       CALL write_trajectory(force_env, root_section, itimes, time*femtoseconds, dt*femtoseconds, md_ener%epot,&
            "FORCES", my_pos, my_act, middle_name="frc", extended_xmol_title=.TRUE., error=error)

       ! Print Simulation Cell 
       CALL write_simulation_cell(cell, motion_section, itimes, time*femtoseconds, my_pos, my_act, error)

       ! Print Thermostats status
       CALL print_thermostats_status(thermostats, para_env, my_pos, my_act, cell, itimes, time, error)

       ! Print Barostat status
       CALL print_barostat_status(barostat, simpar, my_pos, my_act, cell, itimes, time, error)

       ! Print Stress Tensor
       CALL write_stress_tensor(virial, cell, motion_section, unit, itimes, time*femtoseconds,&
            my_pos, my_act, error)

       ! Temperature per Kinds
       IF(simpar%temperature_per_kind) THEN
          tempkind=cp_print_key_unit_nr(logger,motion_section,"MD%PRINT%TEMP_KIND",&
               extension=".temp",file_position=my_pos, file_action=my_act,error=error)
          IF( tempkind > 0 ) THEN
             nkind = SIZE(md_ener%temp_kind)
             fmd="(I10,F20.3,"//TRIM(ADJUSTL(cp_to_string(nkind)))//"F20.9)"     
             fmd=TRIM(fmd)
             WRITE (tempkind,fmd)itimes,time*femtoseconds, md_ener%temp_kind(1:nkind)
             CALL m_flush(tempkind)
          END IF
          CALL cp_print_key_finished_output(tempkind,logger,motion_section,"MD%PRINT%TEMP_KIND", error=error)
       ELSE
          print_key => section_vals_get_subs_vals(motion_section,"MD%PRINT%TEMP_KIND",error=error)
          CALL cp_assert(.NOT.BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file),&
               cp_warning_level,cp_assertion_failed,routineP,&
               "The print_key MD%PRINT%TEMP_KIND has been activated but the "//&
               "calculation of the temperature per kind has not been requested. "//&
               "Please switch on the keyword MD%TEMP_KIND. "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
       END IF

       ! Core/Shell Model
       IF(shell_present) THEN
          CALL force_env_get(force_env, harmonic_shell=harm_shell, error=error)
          CALL cp_subsys_get(subsys, shell_particles=shell_particles, core_particles=core_particles, error=error)

          ! Print Shell Energy
          shene=cp_print_key_unit_nr(logger,motion_section,"MD%PRINT%SHELL_ENERGY",&
                extension=".shener",file_position=my_pos, file_action=my_act, &
                file_form="FORMATTED", is_new_file=new_file,error=error)
          IF (shene>0) THEN
            IF(new_file) THEN
              WRITE (shene,'("#",3X,A,3X,A,3X,3(5X,A,5X))')"Step Nr.","Time[fs]","Kin.[a.u.]",&
                      "Temp.[K]","Pot.[a.u.]" 
            END IF

             WRITE (shene,"(I8,F12.3,F20.9,F20.9,F20.9,F20.9 )")&
                  itimes,time*femtoseconds,md_ener%ekin_shell,md_ener%temp_shell,harm_shell
             CALL m_flush(shene)
          END IF
          CALL cp_print_key_finished_output(shene,logger,motion_section,"MD%PRINT%SHELL_ENERGY", error=error)

          ! Print Shell Coordinates
          CALL write_trajectory (force_env, root_section, itimes, time*femtoseconds, dt*femtoseconds, md_ener%epot,&
               "SHELL_TRAJECTORY", my_pos, my_act, "shpos", shell_particles, extended_xmol_title=.TRUE., error=error)

          IF(shell_adiabatic) THEN
             ! Print Shell Velocities
             CALL write_trajectory (force_env, root_section, itimes, time*femtoseconds, dt*femtoseconds, md_ener%epot,&
                  "SHELL_VELOCITIES", my_pos, my_act, "shvel", shell_particles, extended_xmol_title=.TRUE., error=error)

             ! Print Shell Forces
             CALL write_trajectory (force_env, root_section, itimes, time*femtoseconds, dt*femtoseconds, md_ener%epot,&
                  "SHELL_FORCES", my_pos, my_act, "shfrc", shell_particles, extended_xmol_title=.TRUE., error=error)

             ! Print Core Coordinates
             CALL write_trajectory (force_env, root_section, itimes, time*femtoseconds, dt*femtoseconds, md_ener%epot,&
                  "CORE_TRAJECTORY", my_pos, my_act, "copos", core_particles, extended_xmol_title=.TRUE., error=error)

             ! Print Core Velocities
             CALL write_trajectory (force_env, root_section, itimes, time*femtoseconds, dt*femtoseconds, md_ener%epot,&
                  "CORE_VELOCITIES", my_pos, my_act, "covel", core_particles, extended_xmol_title=.TRUE., error=error)

             ! Print Core Forces
             CALL write_trajectory (force_env, root_section, itimes, time*femtoseconds, dt*femtoseconds, md_ener%epot,&
                  "CORE_FORCES", my_pos, my_act, "cofrc", core_particles, extended_xmol_title=.TRUE., error=error)

             ! Temperature per Kinds
             IF(simpar%temperature_per_kind) THEN
                tempkind=cp_print_key_unit_nr(logger,motion_section,"MD%PRINT%TEMP_SHELL_KIND",&
                     extension=".shtemp", file_position=my_pos, file_action=my_act,error=error)
                IF( tempkind > 0 ) THEN
                   nkind = SIZE(md_ener%temp_shell_kind)
                   fmd="(I10,F20.3,"//TRIM(ADJUSTL(cp_to_string(nkind)))//"F20.9)"     
                   fmd=TRIM(fmd)
                   WRITE (tempkind,fmd)itimes,time*femtoseconds, md_ener%temp_shell_kind(1:nkind)
                   CALL m_flush(tempkind)
                END IF
                CALL cp_print_key_finished_output(tempkind, logger, motion_section,&
                     "MD%PRINT%TEMP_SHELL_KIND", error=error) 
             ELSE
                print_key => section_vals_get_subs_vals(motion_section,"MD%PRINT%TEMP_SHELL_KIND",error=error)
                CALL cp_assert(.NOT.BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file),&
                     cp_warning_level,cp_assertion_failed,routineP,&
                     "The print_key MD%PRINT%TEMP_SHELL_KIND has been activated but the "//&
                     "calculation of the temperature per kind has not been requested. "//&
                     "Please switch on the keyword MD%TEMP_KIND. "//&
CPSourceFileRef,&
                     only_ionode=.TRUE.)
             END IF
          END IF
       END IF
    END IF
    init = .FALSE.
    CALL set_md_env(md_env,init=init,error=error)
    CALL timestop(handle)
  END SUBROUTINE md_write_output

! *****************************************************************************
!> \brief This routine prints all basic information during MD steps
!> \par History
!>      - 10.2008 - Teodoro Laino [tlaino] - University of Zurich
!>                  Refactoring: split into an independent routine.
!>                  All output on screen must be included in this function!
!> \author CJM
! *****************************************************************************
  SUBROUTINE md_write_info_low(simpar, md_ener, unit, qmmm, virial, reftraj, cell,&
       abc, cell_angle, itimes, dt, time, used_time, averages, econs, pv_scalar, &
       pv_xx, hugoniot, nat, init, logger, motion_section, my_pos, my_act, error)

    TYPE(simpar_type), POINTER               :: simpar
    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(unit_convert_type), POINTER         :: unit
    LOGICAL, INTENT(IN)                      :: qmmm
    TYPE(virial_type), POINTER               :: virial
    TYPE(reftraj_type), POINTER              :: reftraj
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: abc, cell_angle
    INTEGER, POINTER                         :: itimes
    REAL(KIND=dp), INTENT(IN)                :: dt
    REAL(KIND=dp), POINTER                   :: time, used_time
    TYPE(average_quantities_type), POINTER   :: averages
    REAL(KIND=dp), INTENT(IN)                :: econs, pv_scalar, pv_xx, &
                                                hugoniot
    INTEGER, INTENT(IN)                      :: nat
    LOGICAL, INTENT(IN)                      :: init
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: motion_section
    CHARACTER(LEN=default_string_length), &
      INTENT(IN)                             :: my_pos, my_act
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'md_write_info_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iw
    LOGICAL                                  :: failure

    failure = .FALSE.
    ! Print to the screen info about MD
    iw = cp_print_key_unit_nr(logger,motion_section,"MD%PRINT%PROGRAM_RUN_INFO",&
         extension=".mdLog",file_position=my_pos,file_action=my_act,error=error)
    ! Performing protocol relevant to the first step of an MD run
    IF (iw>0) THEN
       IF (init) THEN
          ! Write initial values of quantities of interest
          WRITE (iw,*)
          WRITE (iw,'( A )') ' MD_ENERGIES| Initialization proceeding'
          WRITE (iw,*)
          WRITE (iw,'( )')
          WRITE (iw,'( A,A )') ' ******************************** ', &
               'GO CP2K GO! **********************************'
          WRITE (iw,'( A,A,T40,A,T60,1(1X,E20.12) )') ' INITIAL POTENTIAL ENERGY', &
               ADJUSTL( unit%e_label), '= ', md_ener%epot
          IF (simpar%ensemble/=reftraj_ensemble) THEN
             IF (.NOT.qmmm) THEN
                ! NO QM/MM info
                WRITE (iw,'( A,A,T40,A,T60,1(1X,E20.12) )') ' INITIAL KINETIC ENERGY', &
                     ADJUSTL( unit%e_label), '= ', md_ener%ekin
                WRITE (iw,'( A,A,T40,A,T60,1(1X,F20.3) )') ' INITIAL TEMPERATURE', &
                     ADJUSTL( unit%temp_label ), '= ', md_ener%temp_part
             ELSE
                WRITE (iw,'( A,A,T40,A,T60,1(1X,E20.12) )') ' TOTAL INITIAL KINETIC ENERGY', &
                     ADJUSTL( unit%e_label), '= ', md_ener%ekin
                WRITE (iw,'( A,A,T40,A,T60,1(1X,E20.12) )') ' QM INITIAL KINETIC ENERGY', &
                     ADJUSTL( unit%e_label), '= ', md_ener%ekin_qm
                WRITE (iw,'( A,A,T40,A,T60,1(1X,F20.3) )') ' TOTAL INITIAL TEMPERATURE', &
                     ADJUSTL( unit%temp_label ), '= ', md_ener%temp_part
                WRITE (iw,'( A,A,T40,A,T60,1(1X,F20.3) )') ' QM INITIAL TEMPERATURE', &
                     ADJUSTL( unit%temp_label ), '= ', md_ener%temp_qm
             END IF
          END IF
          IF ( simpar%ensemble == nph_uniaxial_ensemble .OR.&
               simpar%ensemble == nph_uniaxial_damped_ensemble.OR.&
               simpar%ensemble == npt_i_ensemble.OR.&
               simpar%ensemble == npt_f_ensemble.OR.&
               simpar%ensemble == npe_i_ensemble.OR.&
               simpar%ensemble == npe_f_ensemble)  &
               WRITE (iw,'( A,A,T40,A,T60,1(1X,E20.12) )') ' INITIAL BAROSTAT TEMP', &
               ADJUSTL( unit%temp_label), '= ', md_ener%temp_baro
          IF ( virial%pv_availability ) &
               WRITE (iw,'( A,A,T40,A,T60,1(1X,E20.12) )') ' INITIAL PRESSURE', &
               ADJUSTL( unit%pres_label), '= ', pv_scalar
          IF ( simpar%ensemble == nph_uniaxial_ensemble .OR.&
               simpar%ensemble == nph_uniaxial_damped_ensemble)&
               WRITE (iw,'( A,A,T40,A,T60,1(1X,E20.12) )') ' INITIAL HUGONIOT CONSTRAINT', &
               ADJUSTL( unit%temp_label), '= ', hugoniot
          IF ( simpar%ensemble == nph_uniaxial_ensemble .OR.&
               simpar%ensemble == nph_uniaxial_damped_ensemble)&
               WRITE (iw,'( A,A,T40,A,T60,1(1X,E20.12) )') ' INITIAL E0', &
               ADJUSTL( unit%e_label), '= ', simpar%e0
          WRITE (iw,'( A,A,T40,A,T60,1(1X,E20.12) )') ' INITIAL VOLUME', &
               ADJUSTL( unit%vol_label), '= ', cell%deth
          WRITE (iw,'( A,A,T29,A,T33,3(1X,E15.7) )') ' INITIAL CELL LNTHS', &
               ADJUSTL(unit%l_label), '= ', abc(1), abc(2), abc(3)
          WRITE (iw,'( A,A,T29,A,T33,3(1X,E15.7) )') ' INITIAL CELL ANGLS', &
               ADJUSTL(unit%angl_label), '= ', cell_angle(3), cell_angle(2), cell_angle(1)
          WRITE (iw,'( A,A )') ' ******************************** ', &
               'GO CP2K GO! **********************************'
       ELSE
          ! Write seuquential values of quantities of interest
          WRITE (iw,'( )')
          WRITE (iw,'( A,A )') ' **************************************', &
               '*****************************************'
          WRITE (iw,'( A,T61,A )') ' ENSEMBLE TYPE                = ', &
               ADJUSTR( ensemble_name_types(simpar%ensemble)(1:20))
          WRITE (iw,'( A,T71,I10 )') ' STEP NUMBER                  = ', itimes
          IF(simpar%variable_dt) THEN
            WRITE (iw,'( A,T60,1(1X,F20.6) )') ' TIME STEP [FS]             = ', dt*femtoseconds
          END IF
          WRITE (iw,'( A,T60,1(1X,F20.6) )') ' TIME [FS]                    = ', time*femtoseconds
          WRITE (iw,'( A,T60,1(1X,E20.12) )') ' CONSERVED QNTY               = ', &
               md_ener%constant
          WRITE (iw,'( )')
          WRITE (iw,'( T42,A )') 'INSTANTANEOUS        AVERAGES'
          WRITE (iw,'( A,T39,2(1X,F20.2) )') ' CPU [S]                      = ', &
               used_time, averages%avecpu
          WRITE (iw,'( A,T39,2(1X,E20.12) )') ' {E-E0}/{k_b*N_at}            = ', &
               econs, averages%econs
          WRITE (iw,'( A,A,T31,A,T39,2(1X,E20.12) )') ' POTENTIAL ENERGY', &
               ADJUSTL(unit%e_label), '= ', md_ener%epot, averages%avepot
          IF(simpar%ensemble/=reftraj_ensemble) THEN
             IF (.NOT.qmmm) THEN
                ! No QM/MM info
                WRITE (iw,'( A,A,T31,A,T39,2(1X,E20.12) )') ' KINETIC ENERGY', &
                     ADJUSTL(unit%e_label), '= ', md_ener%ekin, averages%avekin
                WRITE (iw,'( A,A,T31,A,T39,2(1X,F20.3) )') ' TEMPERATURE', &
                     ADJUSTL(unit%temp_label), '= ', md_ener%temp_part, averages%avetemp
             ELSE
                WRITE (iw,'( A,A,T31,A,T39,2(1X,E20.12) )') ' TOTAL KINETIC ENERGY', &
                     ADJUSTL(unit%e_label), '= ', md_ener%ekin, averages%avekin
                WRITE (iw,'( A,A,T31,A,T39,2(1X,E20.12) )') ' QM KINETIC ENERGY', &
                     ADJUSTL(unit%e_label), '= ', md_ener%ekin_qm, averages%avekin_qm
                WRITE (iw,'( A,A,T31,A,T39,2(1X,F20.3) )') ' TOTAL TEMPERATURE', &
                     ADJUSTL(unit%temp_label), '= ', md_ener%temp_part, averages%avetemp
                WRITE (iw,'( A,A,T31,A,T39,2(1X,F20.3) )') ' QM TEMPERATURE', &
                     ADJUSTL(unit%temp_label), '= ', md_ener%temp_qm, averages%avetemp_qm
             END IF
          END IF
          IF ( virial%pv_availability ) &
               WRITE (iw,'( A,A,T31,A,T39,2(1X,E20.12) )') ' PRESSURE', &
               ADJUSTL(unit%pres_label), '= ', pv_scalar, averages%avepress
          IF ( simpar%ensemble == nph_uniaxial_ensemble .OR.&
               simpar%ensemble == nph_uniaxial_damped_ensemble) THEN
             WRITE (iw,'( A,A,T31,A,T39,2(1X,E20.12) )') ' P_XX', &
                  ADJUSTL(unit%pres_label), '= ', pv_xx, averages%avepxx
             WRITE (iw,'( A,A,T31,A,T39,2(1X,E20.12) )') ' HUGONIOT', &
                  ADJUSTL(unit%temp_label), '= ', hugoniot/3._dp/nat*kelvin,&
                  averages%avehugoniot/3._dp/nat*kelvin
          END IF
          IF ( simpar%ensemble == nph_uniaxial_ensemble .OR.&
               simpar%ensemble == nph_uniaxial_damped_ensemble.OR.&
               simpar%ensemble == npt_i_ensemble.OR.&
               simpar%ensemble == npt_f_ensemble.OR.&
               simpar%ensemble == npe_i_ensemble.OR.&
               simpar%ensemble == npe_f_ensemble) THEN
             WRITE (iw,'( A,A,T31,A,T39,2(1X,E20.12) )') ' BAROSTAT TEMP', &
                  ADJUSTL(unit%temp_label), '= ', md_ener%temp_baro, averages%avetemp_baro
             WRITE (iw,'( A,A,T31,A,T39,2(1X,E20.12) )') ' VOLUME', &
                  ADJUSTL(unit%vol_label), '= ', cell%deth, averages%avevol
             WRITE (iw,'( A,A,T31,A,T33,3(1X,E15.7) )') ' CELL LNTHS', &
                  ADJUSTL(unit%l_label), '= ', abc(1), abc(2), abc(3)
             WRITE (iw,'( A,A,T31,A,T33,3(1X,E15.7) )') ' AVE. CELL LNTHS', &
                  ADJUSTL(unit%l_label), '= ', averages%aveca, averages%avecb, &
                  averages%avecc
          END IF
          IF (simpar%ensemble==npt_f_ensemble .OR. &
              simpar%ensemble == npe_f_ensemble) THEN
             WRITE (iw,'( A,A,T31,A,T33,3(1X,E15.7) )') ' CELL ANGLS', &
                  ADJUSTL(unit%angl_label), '= ', cell_angle(3), cell_angle(2), cell_angle(1)
             WRITE (iw,'( A,A,T31,A,T33,3(1X,E15.7) )') ' AVE. CELL ANGLS', &
                  ADJUSTL(unit%angl_label), '= ', averages%aveal, averages%avebe, averages%avega
          END IF
          IF(simpar%ensemble==reftraj_ensemble) THEN
             IF(reftraj%info%msd) THEN
                IF (reftraj%msd%msd_kind) &
                     WRITE(iw,'(/,A,T50,3f10.5)') ' COM displacement (dx,dy,dz) [angstrom]:  ', &
                     reftraj%msd%drcom(1:3)*angstrom
             END IF
          END IF
          
          WRITE (iw,'( A,A )') ' **************************************', &
               '*****************************************'
          WRITE (iw,'( )')
       END IF
    END IF
    CALL cp_print_key_finished_output(iw,logger,motion_section,&
         "MD%PRINT%PROGRAM_RUN_INFO", error=error)
  END SUBROUTINE md_write_info_low

END MODULE md_energies
