!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!***** cp2k/thermostat_types [1.0] *
!!
!!   NAME
!!     thermostat_types
!!
!!   FUNCTION
!!     Thermostat structure: module containing thermostat available for MD
!!
!!   AUTHOR
!!     teo [tlaino] - University of Zurich - 09.2007
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE thermostat_types
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE csvr_system_types,               ONLY: csvr_system_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_initialize,&
                                             dyn_coeff_set_type
  USE extended_system_types,           ONLY: lnhc_dealloc,&
                                             lnhc_init,&
                                             nhc_set_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE init_extended_system_types,      ONLY: initialize_nhc_baro,&
                                             initialize_nhc_coef,&
                                             initialize_nhc_part,&
                                             initialize_nhc_shell
  USE init_extended_system_variables,  ONLY: init_nhc_forces
  USE input_constants,                 ONLY: do_thermo_nose
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: simulation_parameters_type
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_types_new,              ONLY: global_constraint_type
  USE particle_list_types,             ONLY: particle_list_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: thermostat_type,&
            create_thermostat_type,&
            release_thermostat_type,&
            retain_thermostat_type

  TYPE thermostat_type
     INTEGER :: id_nr, ref_count, type_of_thermostat
     TYPE (nhc_set_type),     POINTER :: nhc
     TYPE (csvr_system_type), POINTER :: csvr
  END TYPE thermostat_type

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'thermostat_types'
  INTEGER, PRIVATE, SAVE :: last_thermostat_id_nr=0,&
                            do_thermo_part  = 0,&
                            do_thermo_kg    = 1,&
                            do_thermo_shell = 2,&
                            do_thermo_baro  = 3

CONTAINS

!!****f* thermostat_types/create_thermostat_type [1.0] *
!!
!!   NAME
!!     create_thermostat_type
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
!!     09.2007 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE create_thermostat_type(thermostat,motion_section, &
       force_env, simpar, para_env, globenv, global_section, error )
    TYPE(thermostat_type), POINTER           :: thermostat
    TYPE(section_vals_type), POINTER         :: motion_section
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: global_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_thermostat_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nshell, stat
    LOGICAL :: check, explicit_baro, explicit_kg, explicit_part, &
      explicit_shell, failure, save_mem, shell_adiabatic, shell_present
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(mol_new_list_type), POINTER         :: molecules_new
    TYPE(particle_list_type), POINTER        :: shell_particles
    TYPE(section_vals_type), POINTER :: kg_coeff_section, print_section, &
      thermo_baro_section, thermo_kg_section, thermo_part_section, &
      thermo_shell_section, work_section

    failure =.FALSE.
    check = .NOT.ASSOCIATED(thermostat)
    CPPrecondition(check,cp_failure_level,routineP,error,failure)
    thermo_part_section  => section_vals_get_subs_vals(motion_section,"MD%THERMOSTAT",error=error)
    thermo_kg_section    => section_vals_get_subs_vals(motion_section,"MD%KG_COEFF%THERMOSTAT",error=error)
    thermo_shell_section => section_vals_get_subs_vals(motion_section,"MD%SHELL%THERMOSTAT",error=error)
    thermo_baro_section  => section_vals_get_subs_vals(motion_section,"MD%BAROSTAT%THERMOSTAT",error=error)

    CALL section_vals_val_get(global_section,"SAVE_MEM",l_val=save_mem,error=error)
    CALL section_vals_get(thermo_part_section, explicit=explicit_part, error=error)
    CALL section_vals_get(thermo_kg_section, explicit=explicit_kg, error=error)
    CALL section_vals_get(thermo_shell_section, explicit=explicit_shell, error=error)
    CALL section_vals_get(thermo_baro_section, explicit=explicit_baro, error=error)

    ALLOCATE(thermostat, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    last_thermostat_id_nr = last_thermostat_id_nr + 1
    thermostat%id_nr = last_thermostat_id_nr
    thermostat%ref_count =  1
    NULLIFY(thermostat%nhc)
    NULLIFY(thermostat%csvr)
    
    thermostat%type_of_thermostat = do_thermo_nose
    ALLOCATE(thermostat%nhc)
    NULLIFY(thermostat%nhc%nhc_part,&
            thermostat%nhc%nhc_baro,&
            thermostat%nhc%nhc_coef,&
            thermostat%nhc%nhc_shell)
    
    ALLOCATE (thermostat%nhc%nhc_part,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (thermostat%nhc%nhc_baro,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (thermostat%nhc%nhc_coef,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (thermostat%nhc%nhc_shell,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL lnhc_init(thermostat%nhc%nhc_part,error=error)
    CALL lnhc_init(thermostat%nhc%nhc_baro,error=error)
    CALL lnhc_init(thermostat%nhc%nhc_coef,error=error)
    CALL lnhc_init(thermostat%nhc%nhc_shell,error=error)

    ! Particles
    CALL force_env_get(force_env, subsys=subsys, cell=cell, error=error )
    CALL cp_subsys_get(subsys, local_molecules_new=local_molecules,&
         molecule_kinds_new=molecule_kinds_new,molecules_new=molecules_new,&
         gci=gci,error=error)
       
    ! Initialize or possibly restart Nose on Particles
    work_section => section_vals_get_subs_vals(thermo_part_section,"NOSE",error=error)
    print_section => section_vals_get_subs_vals(motion_section,"MD%PRINT",error=error)
    CALL initialize_nhc_part(cell, simpar, local_molecules, molecules_new%els,&
         molecule_kinds_new%els, para_env, globenv, thermostat%nhc%nhc_part,&
         print_section=print_section, nose_section=work_section, gci=gci,&
         save_mem=save_mem, error=error)
    CALL init_nhc_forces ( nhc_part=thermostat%nhc%nhc_part,error=error)

    ! Coefficients
    CALL force_env_get(force_env, subsys=subsys, dyn_coeff_set=dyn_coeff_set, error=error )
       
    IF (ASSOCIATED(dyn_coeff_set)) THEN
       ! Initialize coef thermostats and setup coefficient forces, velocities and effective masses
       ! Possibly restart KG coefficients
       kg_coeff_section  => section_vals_get_subs_vals(motion_section,"MD%KG_COEFF",error=error)
       work_section      => section_vals_get_subs_vals(kg_coeff_section,"THERMOSTAT%NOSE",error=error)
       CALL initialize_nhc_coef(dyn_coeff_set, para_env, globenv,&
            thermostat%nhc%nhc_coef, coeff_section=kg_coeff_section,&
            nose_section=work_section, dt=simpar%dt, save_mem=save_mem,&
            error=error)
       
       CALL dyn_coeff_set_initialize(dyn_coeff_set, motion_section, para_env,&
            globenv, coeff_section=kg_coeff_section, error=error)
       CALL init_nhc_forces (nhc_coef=thermostat%nhc%nhc_coef,error=error)
    ELSE
       DEALLOCATE (thermostat%nhc%nhc_coef,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    
    ! Core-Shell Model
    CALL force_env_get(force_env, subsys=subsys, error=error )
    
    CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds,&
         local_particles=local_particles,shell_particles=shell_particles,&
         error=error)
    
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kinds%els,&
         shell_present=shell_present,shell_adiabatic=shell_adiabatic)
    
    IF(shell_present .AND. shell_adiabatic) THEN
       ! Initialize SHELL thermostats, possibly restart
       nshell = SIZE(shell_particles%els)
       work_section => section_vals_get_subs_vals(thermo_shell_section,"NOSE",error=error)
       CALL initialize_nhc_shell(simpar, para_env, globenv, thermostat%nhc%nhc_shell,&
            atomic_kinds%els, local_particles, work_section, save_mem,&
            error=error)
       CALL init_nhc_forces( nhc_shell=thermostat%nhc%nhc_shell,error=error)
    ELSE IF(simpar%shell_nose .AND. .NOT.  shell_adiabatic) THEN                            
       CALL stop_program(routineN,moduleN,__LINE__,&
            "Nose on the core-shell motion only with adiabatic shell-model ") 
    ELSE
       DEALLOCATE (thermostat%nhc%nhc_shell,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    
    ! Barostat
    CALL force_env_get( force_env, cell=cell, error=error )
    ! Initialize and possibly restart thermostat of barostat
    work_section => section_vals_get_subs_vals(thermo_baro_section,"NOSE",error=error)
    CALL initialize_nhc_baro(simpar, para_env, globenv, thermostat%nhc%nhc_baro,&
         nose_section=work_section, save_mem=save_mem, error=error)
    CALL init_nhc_forces ( nhc_baro=thermostat%nhc%nhc_baro,error=error)
    
  END SUBROUTINE create_thermostat_type

!!****f* thermostat_types/retain_thermostat_type [1.0] *
!!
!!   NAME
!!     retain_thermostat_type
!!
!!   FUNCTION
!!     retains the given thermostat
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     -
!!
!!   AUTHOR
!!    Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [tlaino]
!! 
!!*** **********************************************************************
  SUBROUTINE retain_thermostat_type(thermostat,error)
    TYPE(thermostat_type), POINTER           :: thermostat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'retain_thermostat_type', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    
    IF (.NOT. failure) THEN
       IF (ASSOCIATED(thermostat)) THEN
          CPPrecondition(thermostat%ref_count>0,cp_failure_level,routineP,error,failure)
          thermostat%ref_count=thermostat%ref_count+1
       END IF
    END IF

  END SUBROUTINE retain_thermostat_type

!!****f* thermostat_types/release_thermostat_type [1.0] *
!!
!!   NAME
!!     release_thermostat_type
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
!!     09.2007 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE release_thermostat_type(thermostat, error)
    TYPE(thermostat_type), POINTER           :: thermostat
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_thermostat_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: check, failure

    failure =.FALSE.
    check = ASSOCIATED(thermostat)
    IF (check) THEN
       check = thermostat%ref_count>0
       CPPrecondition(check,cp_failure_level,routineP,error,failure)
       thermostat%ref_count=thermostat%ref_count-1
       IF (thermostat%ref_count<1) THEN
          IF (ASSOCIATED(thermostat%nhc)) THEN
             CALL lnhc_dealloc(thermostat%nhc%nhc_part,error=error)
             CALL lnhc_dealloc(thermostat%nhc%nhc_baro,error=error)
             CALL lnhc_dealloc(thermostat%nhc%nhc_coef,error=error)
             CALL lnhc_dealloc(thermostat%nhc%nhc_shell,error=error)
             DEALLOCATE(thermostat%nhc, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(thermostat%csvr)) THEN
             DEALLOCATE(thermostat%csvr, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          DEALLOCATE(thermostat, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE release_thermostat_type
!!*****
!******************************************************************************

END MODULE thermostat_types
!******************************************************************************
