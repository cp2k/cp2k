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
  USE atomic_kind_types,               ONLY: get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE init_extended_system_types,      ONLY: initialize_nhc_baro,&
                                             initialize_nhc_coef,&
                                             initialize_nhc_part,&
                                             initialize_nhc_shell
  USE init_extended_system_variables,  ONLY: init_nhc_forces
  USE input_constants,                 ONLY: npt_f_ensemble,&
                                             npt_i_ensemble
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: simulation_parameters_type
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_types_new,              ONLY: global_constraint_type
  USE termination,                     ONLY: stop_program
  USE thermostat_types,                ONLY: allocate_thermostats,&
                                             create_thermostat_type,&
                                             release_thermostat_info,&
                                             release_thermostat_type,&
                                             release_thermostats,&
                                             thermostats_type
  USE thermostat_utils,                ONLY: compute_degrees_of_freedom,&
                                             setup_thermostat_info
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: create_thermostats

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

    LOGICAL :: apply_thermo_baro, do_thermo_shell, explicit_baro, &
      explicit_coef, explicit_part, explicit_shell, failure, save_mem, &
      shell_adiabatic, shell_present
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

    ! Compute Degrees of Freedom
    CALL cp_subsys_get(subsys, molecule_kinds_new=molecule_kinds_new, local_molecules_new=local_molecules,&
            gci=gci, error=error) 
    CALL compute_degrees_of_freedom(thermostats, cell, simpar, molecule_kinds_new%els,&
         local_molecules, print_section, gci, error)

    ! Particles
    IF (explicit_part) THEN
       CALL create_thermostat_type(thermostats%thermostat_part, thermo_part_section, error)
       CALL cp_subsys_get(subsys, molecules_new=molecules_new, error=error)
       
       ! Initialize or possibly restart Nose on Particles
       work_section => section_vals_get_subs_vals(thermo_part_section,"NOSE",error=error)
       CALL initialize_nhc_part(thermostats%thermostat_info_part, simpar, local_molecules,&
            molecules_new%els, molecule_kinds_new%els, para_env, globenv,&
            thermostats%thermostat_part%nhc, nose_section=work_section,&
            gci=gci, save_mem=save_mem, error=error)
       CALL init_nhc_forces ( nhc=thermostats%thermostat_part%nhc,error=error)
    END IF

    ! Core-Shell Model
    IF (explicit_shell) THEN
       CALL create_thermostat_type(thermostats%thermostat_shell, thermo_shell_section, error)
       CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds, error=error)
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kinds%els,shell_present=shell_present,&
            shell_adiabatic=shell_adiabatic)
       IF(simpar%shell_nose .AND. .NOT.  shell_adiabatic) THEN                            
          CALL stop_program(routineN,moduleN,__LINE__,&
               "Nose on the core-shell motion only with adiabatic shell-model ")
       END IF
       do_thermo_shell = (shell_present .AND. shell_adiabatic)
       CALL setup_thermostat_info(thermostats%thermostat_info_shell, molecule_kinds_new%els,&
            local_molecules, simpar%shell_region, simpar%ensemble, shell=do_thermo_shell,&
            error=error)
       IF(do_thermo_shell) THEN
          ! Initialize SHELL thermostats, possibly restart
          work_section => section_vals_get_subs_vals(thermo_shell_section,"NOSE",error=error)
          CALL initialize_nhc_shell(thermostats%thermostat_info_shell, simpar, local_molecules,&
            molecules_new%els, molecule_kinds_new%els, para_env, globenv,&
            thermostats%thermostat_shell%nhc, nose_section=work_section,gci=gci,&
            save_mem=save_mem, error=error) 
          CALL init_nhc_forces( nhc=thermostats%thermostat_shell%nhc,error=error)
       ELSE
          CALL release_thermostat_type(thermostats%thermostat_shell, error)
          CALL release_thermostat_info(thermostats%thermostat_info_shell, error)
       END IF
    END IF

    ! Coefficients
    IF (explicit_coef) THEN
       CALL create_thermostat_type(thermostats%thermostat_coef, thermo_coef_section, error)
       CALL force_env_get(force_env, dyn_coeff_set=dyn_coeff_set, error=error )
       IF (ASSOCIATED(dyn_coeff_set)) THEN
          kg_coeff_section  => section_vals_get_subs_vals(motion_section,"MD%KG_COEFF",error=error)
          work_section      => section_vals_get_subs_vals(kg_coeff_section,"THERMOSTAT%NOSE",error=error)
          CALL initialize_nhc_coef(dyn_coeff_set, para_env, globenv,&
               thermostats%thermostat_coef%nhc, coeff_section=kg_coeff_section,&
               nose_section=work_section, dt=simpar%dt, save_mem=save_mem,&
               error=error)
          CALL init_nhc_forces (nhc=thermostats%thermostat_coef%nhc,error=error)
       ELSE
          CALL release_thermostat_type(thermostats%thermostat_coef, error)
       END IF
    END IF

    ! Barostat
    IF (apply_thermo_baro) THEN
       CALL create_thermostat_type(thermostats%thermostat_baro, thermo_baro_section, error)
       work_section => section_vals_get_subs_vals(thermo_baro_section,"NOSE",error=error)
       ! Initialize and possibly restart thermostat of barostat
       CALL initialize_nhc_baro(simpar, para_env, globenv, thermostats%thermostat_baro%nhc,&
            nose_section=work_section, save_mem=save_mem, error=error)
       CALL init_nhc_forces ( nhc=thermostats%thermostat_baro%nhc, error=error)
    END IF

    ! If no thermostats have been allocated deallocate the full structure
    IF(  (.NOT.ASSOCIATED(thermostats%thermostat_part)).AND.&
         (.NOT.ASSOCIATED(thermostats%thermostat_shell)).AND.&
         (.NOT.ASSOCIATED(thermostats%thermostat_coef)).AND.&
         (.NOT.ASSOCIATED(thermostats%thermostat_baro)) ) THEN
       CALL release_thermostats(thermostats, error)
    END IF

  END SUBROUTINE create_thermostats

END MODULE thermostat_methods
!******************************************************************************
