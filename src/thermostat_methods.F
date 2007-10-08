!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!***** cp2k/thermostat_methods [1.0] *
!!
!!   NAME
!!     thermostat_methods
!!
!!   FUNCTION
!!     Methods for Thermostats 
!!
!!   AUTHOR
!!     teo [tlaino] - University of Zurich - 10.2007
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE thermostat_methods
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_type
  USE extended_system_dynamics,        ONLY: lnhc_barostat,&
                                             lnhc_coefs,&
                                             lnhc_particle_set,&
                                             lnhc_particles,&
                                             lnhc_shell_set,&
                                             lnhc_shells
  USE extended_system_init,            ONLY: initialize_nhc_baro,&
                                             initialize_nhc_coef,&
                                             initialize_nhc_part,&
                                             initialize_nhc_shell
  USE extended_system_types,           ONLY: npt_info_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: do_thermo_csvr,&
                                             do_thermo_nose,&
                                             do_thermo_same_as_part,&
                                             npt_f_ensemble,&
                                             npt_i_ensemble,&
                                             nve_ensemble,&
                                             nvt_ensemble
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: simulation_parameters_type
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: global_constraint_type,&
                                             molecule_type
  USE particle_types,                  ONLY: particle_type
  USE thermostat_types,                ONLY: allocate_thermostats,&
                                             create_thermostat_type,&
                                             release_thermostat_info,&
                                             release_thermostat_type,&
                                             release_thermostats,&
                                             thermostat_type,&
                                             thermostats_type
  USE thermostat_utils,                ONLY: compute_degrees_of_freedom,&
                                             setup_thermostat_info
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: create_thermostats,&
            apply_thermostat_coeff,&
            apply_thermostat_baro,&
            apply_thermostat_part_begin,&
            apply_thermostat_part_end,&
            apply_thermostat_shell_begin,&
            apply_thermostat_shell_end

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'thermostat_methods'

CONTAINS

!!****f* thermostat_methods/create_thermostats [1.0] *
!!
!!   NAME
!!     create_thermostats
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
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2007 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE create_thermostats(thermostats,motion_section, &
       force_env, simpar, para_env, globenv, global_section, error )
    TYPE(thermostats_type), POINTER          :: thermostats
    TYPE(section_vals_type), POINTER         :: motion_section
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: global_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_thermostats', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: thermostat_type
    LOGICAL :: apply_general_thermo, apply_thermo_baro, apply_thermo_coeff, &
      do_thermo_shell, explicit_baro, explicit_coef, explicit_part, &
      explicit_shell, failure, save_mem, shell_adiabatic, shell_present
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(mol_new_list_type), POINTER         :: molecules_new
    TYPE(section_vals_type), POINTER :: kg_coeff_section, print_section, &
      thermo_baro_section, thermo_coef_section, thermo_part_section, &
      thermo_shell_section, work_section

    failure =.FALSE.
    CALL allocate_thermostats(thermostats, error)
    thermo_part_section  => section_vals_get_subs_vals(motion_section,"MD%THERMOSTAT",error=error)
    thermo_coef_section  => section_vals_get_subs_vals(motion_section,"MD%KG_COEFF%THERMOSTAT",error=error)
    thermo_shell_section => section_vals_get_subs_vals(motion_section,"MD%SHELL%THERMOSTAT",error=error)
    thermo_baro_section  => section_vals_get_subs_vals(motion_section,"MD%BAROSTAT%THERMOSTAT",error=error)
    print_section => section_vals_get_subs_vals(motion_section,"MD%PRINT",error=error)

    CALL force_env_get(force_env, subsys=subsys, cell=cell, error=error )
    CALL section_vals_val_get(global_section,"SAVE_MEM",l_val=save_mem,  error=error)
    CALL section_vals_get(thermo_part_section,  explicit=explicit_part,  error=error)
    CALL section_vals_get(thermo_coef_section,  explicit=explicit_coef,    error=error)
    CALL section_vals_get(thermo_shell_section, explicit=explicit_shell, error=error)
    CALL section_vals_get(thermo_baro_section,  explicit=explicit_baro,  error=error)

    apply_thermo_baro = (simpar%ensemble == npt_f_ensemble).OR.&
                        (simpar%ensemble == npt_i_ensemble)

    apply_thermo_coeff= (simpar%ensemble == nve_ensemble).OR.&
                        (simpar%ensemble == nvt_ensemble)

    apply_general_thermo = apply_thermo_baro .OR. (simpar%ensemble == nvt_ensemble)

    ! Compute Degrees of Freedom
    CALL cp_subsys_get(subsys, molecule_kinds_new=molecule_kinds_new, local_molecules_new=local_molecules,&
            gci=gci, error=error) 
    CALL compute_degrees_of_freedom(thermostats, cell, simpar, molecule_kinds_new%els,&
         local_molecules, print_section, gci, error)

    ! Particles
    IF (explicit_part) THEN
       IF (apply_general_thermo) THEN
          CALL create_thermostat_type(thermostats%thermostat_part, thermo_part_section, error)
          CALL cp_subsys_get(subsys, molecules_new=molecules_new, error=error)
          ! Initialize thermostat
          IF      (thermostats%thermostat_part%type_of_thermostat==do_thermo_nose) THEN
             ! Initialize or possibly restart Nose on Particles
             work_section => section_vals_get_subs_vals(thermo_part_section,"NOSE",error=error)
             CALL initialize_nhc_part(thermostats%thermostat_info_part, simpar, local_molecules,&
                  molecules_new%els, molecule_kinds_new%els, para_env, globenv,&
                  thermostats%thermostat_part%nhc, nose_section=work_section,&
                  gci=gci, save_mem=save_mem, error=error)
          ELSE IF (thermostats%thermostat_part%type_of_thermostat==do_thermo_csvr) THEN
             ! Initialize or possibly restart CSVR thermostat on Particles
             work_section => section_vals_get_subs_vals(thermo_part_section,"CSVR",error=error)
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END IF
       ELSE
          CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
               "Thermostat for Particles has been defined but the ensemble provided "//&
               "does not support thermostat for Particles! Ignoring thermostat input."//&
CPSourceFileRef,&
               only_ionode=.TRUE.)          
       END IF
    END IF

    ! Core-Shell Model
    IF (explicit_shell) THEN
       IF (apply_general_thermo) THEN
          CALL create_thermostat_type(thermostats%thermostat_shell, thermo_shell_section, error)
          CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds, error=error)
          CALL get_atomic_kind_set(atomic_kind_set=atomic_kinds%els,shell_present=shell_present,&
               shell_adiabatic=shell_adiabatic)
          do_thermo_shell = (shell_present .AND. shell_adiabatic)
          CALL setup_thermostat_info(thermostats%thermostat_info_shell, molecule_kinds_new%els,&
               local_molecules, simpar%shell_region, simpar%ensemble, shell=do_thermo_shell,&
               error=error)
          IF(do_thermo_shell) THEN
             ! Initialize thermostat
             IF      (thermostats%thermostat_shell%type_of_thermostat==do_thermo_nose) THEN
                ! Initialize or possibly restart Nose on Shells
                work_section => section_vals_get_subs_vals(thermo_shell_section,"NOSE",error=error)
                CALL initialize_nhc_shell(thermostats%thermostat_info_shell, simpar, local_molecules,&
                     molecules_new%els, molecule_kinds_new%els, para_env, globenv,&
                     thermostats%thermostat_shell%nhc, nose_section=work_section,gci=gci,&
                     save_mem=save_mem, error=error) 
             ELSE IF (thermostats%thermostat_shell%type_of_thermostat==do_thermo_csvr) THEN
                ! Initialize or possibly restart CSVR thermostat on Shells
                work_section => section_vals_get_subs_vals(thermo_shell_section,"CSVR",error=error)
                CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END IF
          ELSE
             IF (.NOT.shell_adiabatic) THEN
                CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                     "Thermostat for Core-Shell motion only with adiabatic shell-model. "//&
                     "Continuing calculation ignoring the thermostat info! No Thermostat "//&
                     "applied to Shells!"//&
CPSourceFileRef,&
                     only_ionode=.TRUE.) 
             ELSE
                CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                     "Thermostat for Shells has been defined but the system provided "//&
                     "does not contain any Shells! Ignoring thermostat input."//&
CPSourceFileRef,&
                     only_ionode=.TRUE.) 
             END IF
             CALL release_thermostat_type(thermostats%thermostat_shell, error)
             CALL release_thermostat_info(thermostats%thermostat_info_shell, error)
          END IF
       ELSE
          CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
               "Thermostat for Shells has been defined but the ensemble provided "//&
               "does not support thermostat for Shells! Ignoring thermostat input."//&
CPSourceFileRef,&
               only_ionode=.TRUE.)          
       END IF
    END IF

    ! Coefficients
    IF (explicit_coef) THEN
       IF (apply_thermo_coeff) THEN
          CALL create_thermostat_type(thermostats%thermostat_coef, thermo_coef_section, error)
          CALL force_env_get(force_env, dyn_coeff_set=dyn_coeff_set, error=error )
          IF (ASSOCIATED(dyn_coeff_set)) THEN
             kg_coeff_section  => section_vals_get_subs_vals(motion_section,"MD%KG_COEFF",error=error)
             ! Initialize thermostat
             IF      (thermostats%thermostat_coef%type_of_thermostat==do_thermo_nose) THEN
                ! Initialize or possibly restart Nose on Coeffs
                work_section      => section_vals_get_subs_vals(kg_coeff_section,"THERMOSTAT%NOSE",error=error)
                CALL initialize_nhc_coef(dyn_coeff_set, para_env, globenv,&
                     thermostats%thermostat_coef%nhc, coeff_section=kg_coeff_section,&
                     nose_section=work_section, dt=simpar%dt, save_mem=save_mem,&
                     error=error)
             ELSE IF (thermostats%thermostat_coef%type_of_thermostat==do_thermo_csvr) THEN
                ! Initialize or possibly restart CSVR thermostat on Coeffs
                work_section      => section_vals_get_subs_vals(kg_coeff_section,"THERMOSTAT%CSVR",error=error)
                CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END IF
          ELSE
             CALL release_thermostat_type(thermostats%thermostat_coef, error)
          END IF
       ELSE
          CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
               "Thermostat for KG Coefficients has been defined but the ensemble provided "//&
               "does not support thermostat for KG coefficients! Ignoring thermostat input."//&
CPSourceFileRef,&
               only_ionode=.TRUE.)          
       END IF
    END IF

    ! Barostat
    IF (apply_thermo_baro) THEN
       ! Check if we use the same thermostat as particles
       CALL section_vals_val_get(thermo_baro_section,"TYPE",i_val=thermostat_type,error=error)
       work_section => thermo_baro_section
       IF (thermostat_type==do_thermo_same_as_part) work_section => thermo_part_section
       CALL create_thermostat_type(thermostats%thermostat_baro, work_section, error)
       ! Initialize thermostat
       IF      (thermostats%thermostat_baro%type_of_thermostat==do_thermo_nose) THEN
          ! Initialize or possibly restart Nose on Barostat
          work_section => section_vals_get_subs_vals(thermo_baro_section,"NOSE",error=error)
          CALL initialize_nhc_baro(simpar, para_env, globenv, thermostats%thermostat_baro%nhc,&
               nose_section=work_section, save_mem=save_mem, error=error)
       ELSE IF (thermostats%thermostat_baro%type_of_thermostat==do_thermo_csvr) THEN
          ! Initialize or possibly restart CSVR thermostat on Barostat
          work_section => section_vals_get_subs_vals(thermo_baro_section,"CSVR",error=error)
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
    ELSE
       IF (explicit_baro) THEN
          CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
               "Thermostat for Barostat has been defined but the ensemble provided "//&
               "does not support thermostat for Barostat! Ignoring thermostat input."//&
CPSourceFileRef,&
               only_ionode=.TRUE.)          
       END IF
    END IF

    ! If no thermostats have been allocated deallocate the full structure
    IF(  (.NOT.ASSOCIATED(thermostats%thermostat_part)).AND.&
         (.NOT.ASSOCIATED(thermostats%thermostat_shell)).AND.&
         (.NOT.ASSOCIATED(thermostats%thermostat_coef)).AND.&
         (.NOT.ASSOCIATED(thermostats%thermostat_baro)) ) THEN
       CALL release_thermostats(thermostats, error)
    END IF

  END SUBROUTINE create_thermostats

!!****f* thermostat_methods/apply_thermostat_coeff [1.0] *
!!
!!   NAME
!!     apply_thermostat_coeff
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
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2007 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE apply_thermostat_coeff(thermostat, dyn_coeff_set, group, error)
    TYPE(thermostat_type), POINTER           :: thermostat
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    INTEGER, INTENT(IN)                      :: group
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_thermostat_coeff', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(thermostat)) THEN
       IF      (thermostat%type_of_thermostat==do_thermo_nose) THEN
          ! Apply Nose-Hoover Thermostat
          CPPostcondition(ASSOCIATED(thermostat%nhc),cp_failure_level,routineP,error,failure)
          CALL lnhc_coefs(thermostat%nhc, dyn_coeff_set, group)
       ELSE IF (thermostat%type_of_thermostat==do_thermo_csvr) THEN  
          ! Apply CSVR Thermostat
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE apply_thermostat_coeff

!!****f* thermostat_methods/apply_thermostat_baro [1.0] *
!!
!!   NAME
!!     apply_thermostat_baro
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
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2007 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE apply_thermostat_baro(thermostat, npt, group, error)
    TYPE(thermostat_type), POINTER           :: thermostat
    TYPE(npt_info_type), DIMENSION(:, :), &
      POINTER                                :: npt
    INTEGER, INTENT(IN)                      :: group
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_thermostat_baro', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(thermostat)) THEN
       IF      (thermostat%type_of_thermostat==do_thermo_nose) THEN
          ! Apply Nose-Hoover Thermostat
          CPPostcondition(ASSOCIATED(thermostat%nhc),cp_failure_level,routineP,error,failure)
          CALL lnhc_barostat ( thermostat%nhc, npt, group )
       ELSE IF (thermostat%type_of_thermostat==do_thermo_csvr) THEN  
          ! Apply CSVR Thermostat
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE apply_thermostat_baro

!!****f* thermostat_methods/apply_thermostat_part_begin [1.0] *
!!
!!   NAME
!!     apply_thermostat_part_begin
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
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2007 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE apply_thermostat_part_begin(thermostat, molecule_kind_set, molecule_set,&
       particle_set, local_molecules, group, shell_adiabatic, shell_particle_set,&
       core_particle_set, error)

    TYPE(thermostat_type), POINTER           :: thermostat
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    INTEGER, INTENT(IN)                      :: group
    LOGICAL, INTENT(IN), OPTIONAL            :: shell_adiabatic
    TYPE(particle_type), OPTIONAL, POINTER   :: shell_particle_set( : ), &
                                                core_particle_set( : )
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_thermostat_part_begin', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(thermostat)) THEN
       IF      (thermostat%type_of_thermostat==do_thermo_nose) THEN
          ! Apply Nose-Hoover Thermostat
          CPPostcondition(ASSOCIATED(thermostat%nhc),cp_failure_level,routineP,error,failure)
          CALL lnhc_particle_set (thermostat%nhc,molecule_kind_set,molecule_set,&
               particle_set,local_molecules,group,shell_adiabatic,&
               shell_particle_set,core_particle_set)
       ELSE IF (thermostat%type_of_thermostat==do_thermo_csvr) THEN  
          ! Apply CSVR Thermostat
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE apply_thermostat_part_begin

!!****f* thermostat_methods/apply_thermostat_part_end [1.0] *
!!
!!   NAME
!!     apply_thermostat_part_end
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
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2007 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE apply_thermostat_part_end(thermostat, molecule_kind_set, molecule_set,&
       particle_set, local_molecules, vel, group, shell_adiabatic, shell_vel, core_vel, error )

    TYPE(thermostat_type), POINTER           :: thermostat
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    INTEGER, INTENT(IN)                      :: group
    LOGICAL, INTENT(IN), OPTIONAL            :: shell_adiabatic
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL   :: shell_vel(:,:), core_vel(:,:)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_thermostat_part_end', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(thermostat)) THEN
       IF      (thermostat%type_of_thermostat==do_thermo_nose) THEN
          ! Apply Nose-Hoover Thermostat
          CPPostcondition(ASSOCIATED(thermostat%nhc),cp_failure_level,routineP,error,failure)
          CALL lnhc_particles(thermostat%nhc, molecule_kind_set,molecule_set,&
               particle_set, local_molecules, vel, group,  shell_adiabatic, shell_vel, &
               core_vel)
       ELSE IF (thermostat%type_of_thermostat==do_thermo_csvr) THEN  
          ! Apply CSVR Thermostat
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE apply_thermostat_part_end

!!****f* thermostat_methods/apply_thermostat_shell_begin [1.0] *
!!
!!   NAME
!!     apply_thermostat_shell_begin
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
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2007 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE apply_thermostat_shell_begin(thermostat,atomic_kind_set, particle_set,&
       local_particles, shell_particle_set, core_particle_set, group, error)

    TYPE(thermostat_type), POINTER           :: thermostat
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), POINTER             :: shell_particle_set(:), &
                                                core_particle_set(:)
    INTEGER, INTENT(IN)                      :: group
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_thermostat_shell_begin', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(thermostat)) THEN
       IF      (thermostat%type_of_thermostat==do_thermo_nose) THEN
          ! Apply Nose-Hoover Thermostat
          CPPostcondition(ASSOCIATED(thermostat%nhc),cp_failure_level,routineP,error,failure)
          CALL lnhc_shell_set(thermostat%nhc, atomic_kind_set, particle_set, local_particles, &
               shell_particle_set, core_particle_set, group)

       ELSE IF (thermostat%type_of_thermostat==do_thermo_csvr) THEN  
          ! Apply CSVR Thermostat
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE apply_thermostat_shell_begin

!!****f* thermostat_methods/apply_thermostat_shell_end [1.0] *
!!
!!   NAME
!!     apply_thermostat_shell_end
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
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2007 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE apply_thermostat_shell_end(thermostat, atomic_kind_set, particle_set,&
       local_particles, vel, shell_vel, core_vel, group, error)

    TYPE(thermostat_type), POINTER           :: thermostat
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    REAL(KIND=dp), INTENT(INOUT)             :: vel(:,:), shell_vel(:,:), &
                                                core_vel(:,:)
    INTEGER, INTENT(IN)                      :: group
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_thermostat_shell_end', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(thermostat)) THEN
       IF      (thermostat%type_of_thermostat==do_thermo_nose) THEN
          ! Apply Nose-Hoover Thermostat
          CPPostcondition(ASSOCIATED(thermostat%nhc),cp_failure_level,routineP,error,failure)
          CALL lnhc_shells(thermostat%nhc, atomic_kind_set, particle_set, local_particles, &
               vel, shell_vel, core_vel, group)

       ELSE IF (thermostat%type_of_thermostat==do_thermo_csvr) THEN  
          ! Apply CSVR Thermostat
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE apply_thermostat_shell_end

END MODULE thermostat_methods
!******************************************************************************
