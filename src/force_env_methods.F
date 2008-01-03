!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/force_env_types *
!!
!!   NAME
!!     force_env_types
!!
!!   FUNCTION
!!     Interface for the force calculations
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     cjm, FEB-20-2001: pass variable box_ref
!!     cjm, SEPT-12-2002: major reorganization
!!     fawzi, APR-12-2003: introduced force_env (based on the work by CJM&JGH)
!!     fawzi, NOV-3-2004: reorganized interface for f77 interface
!!
!!   SOURCE
!******************************************************************************

MODULE force_env_methods
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE cell_types,                      ONLY: cell_type,&
                                             compare_cells
  USE constraint,                      ONLY: rattle_control,&
                                             shake_control
  USE constraint_fxd,                  ONLY: fix_atom_control
  USE constraint_util,                 ONLY: getold
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_env,                     ONLY: cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_result_types,                 ONLY: cp_result_create,&
                                             cp_result_type,&
                                             erase_result,&
                                             transfer_result_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_p_type,&
                                             cp_subsys_retain,&
                                             cp_subsys_type
  USE cpot_types,                      ONLY: cpot_calc
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_retain
  USE eip_environment_types,           ONLY: eip_env_get,&
                                             eip_env_retain,&
                                             eip_environment_type
  USE eip_silicon,                     ONLY: eip_bazant,&
                                             eip_lenosky
  USE ep_types,                        ONLY: ep_env_calc_e_f,&
                                             ep_env_create,&
                                             ep_env_release,&
                                             ep_env_retain,&
                                             ep_env_type
  USE ewald_environment_types,         ONLY: ewald_environment_type
  USE ewald_pw_methods,                ONLY: ewald_pw_grid_change
  USE ewald_pw_types,                  ONLY: ewald_pw_type
  USE f77_blas
  USE fist_environment_types,          ONLY: fist_env_retain,&
                                             fist_environment_type,&
                                             get_fist_env,&
                                             set_fist_env
  USE fist_force,                      ONLY: fist_force_control
  USE force_env_types,                 ONLY: &
       force_env_get, force_env_get_natom, force_env_p_type, &
       force_env_pos_get, force_env_type, use_eip_force, use_ep_force, &
       use_fist_force, use_kg_force, use_mixed_force, use_prog_name, &
       use_qmmm, use_qs_force
  USE fp_methods,                      ONLY: fp_eval
  USE fparser,                         ONLY: EvalErrType,&
                                             evalf,&
                                             evalfd,&
                                             finalizef,&
                                             initf,&
                                             parsef
  USE global_types,                    ONLY: global_environment_type,&
                                             globenv_retain
  USE input_constants,                 ONLY: &
       do_stress_analytical, do_stress_diagonal_anal, &
       do_stress_diagonal_numer, do_stress_none, do_stress_numerical, &
       dump_xmol, mix_coupled, mix_generic, mix_linear_combination, &
       mix_minimum, mix_restrained, qs_debug_run, use_bazant_eip, &
       use_lenosky_eip
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_retain,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kg_energy,                       ONLY: kg_energies
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_env_retain,&
                                             kg_environment_type,&
                                             set_kg_env
  USE kg_force,                        ONLY: kg_force_control
  USE kg_rho_methods,                  ONLY: kg_rho_create
  USE kg_rho_types,                    ONLY: kg_rho_type
  USE kg_rspw_methods,                 ONLY: kg_rspw_rebuild
  USE kg_rspw_types,                   ONLY: kg_rspw_type
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE message_passing,                 ONLY: mp_sum,&
                                             mp_sync
  USE metadynamics_types,              ONLY: meta_env_retain,&
                                             meta_env_type
  USE mixed_energy_types,              ONLY: mixed_energy_type,&
                                             mixed_force_type
  USE mixed_environment_types,         ONLY: get_mixed_env,&
                                             mixed_env_retain,&
                                             mixed_environment_type
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_types_new,              ONLY: global_constraint_type
  USE particle_list_types,             ONLY: particle_list_p_type,&
                                             particle_list_type
  USE particle_types,                  ONLY: particle_type,&
                                             update_particle_set,&
                                             write_particle_coordinates
  USE qmmm_gpw_energy,                 ONLY: qmmm_el_coupling
  USE qmmm_gpw_forces,                 ONLY: qmmm_forces
  USE qmmm_links_methods,              ONLY: qmmm_added_chrg_coord,&
                                             qmmm_added_chrg_forces,&
                                             qmmm_link_Imomm_coord,&
                                             qmmm_link_Imomm_forces
  USE qmmm_types,                      ONLY: fist_subsys,&
                                             qmmm_env_qm_retain,&
                                             qmmm_env_qm_type,&
                                             qmmm_links_type,&
                                             qs_subsys
  USE qmmm_util,                       ONLY: apply_qmmm_translate
  USE qs_energy,                       ONLY: qs_energies
  USE qs_environment_methods,          ONLY: qs_env_rebuild_pw_env
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_env_retain,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force,                        ONLY: qs_forces
  USE qs_ks_methods,                   ONLY: qs_ks_create
  USE qs_ks_qmmm_methods,              ONLY: ks_qmmm_env_rebuild
  USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                             qs_ks_release
  USE qs_numerical_pressure,           ONLY: qs_calc_numerical_pressure
  USE restraint,                       ONLY: restraint_control
  USE string_utilities,                ONLY: compress
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: cp_virial,&
                                             virial_create,&
                                             virial_p_type,&
                                             virial_release,&
                                             virial_retain,&
                                             virial_set,&
                                             virial_type,&
                                             zero_virial
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'force_env_methods'

  PUBLIC :: force_env_create,&
            ep_create_force_env,&
            force_env_calc_energy_force, &
            force_env_set_cell,&
            force_env_shake,&
            force_env_rattle

  INTEGER, SAVE, PRIVATE :: last_force_env_id=0

!!***
!******************************************************************************

CONTAINS

!!****** force_env_types/force_env_calc_energy_force [1.0] *
!!
!!   NAME
!!     force_env_calc_energy_force
!!
!!   FUNCTION
!!     Interface routine for force and energy calculations
!!
!!   ARGUMENTS
!!     - force_env: the force_env of which you want the energy and forces
!!     - calc_force: if false the forces *might* be left unchanged
!!       or be unvalid, no guarantee on them is done.Defaults to true
!!     - consistent_energies: Performs an additional qs_ks_update_qs_env, so
!!         that the energies are appropriate to the forces, they are in the
!!         non-selfconsistent case not consistent to each other! [08.2005, TdK]
!!
!!   AUTHOR
!!     CJM & fawzi
!!
!!   SOURCE
!******************************************************************************
  RECURSIVE SUBROUTINE force_env_calc_energy_force ( force_env, calc_force, &
       consistent_energies, skip_external_control, eval_energy_forces, error)

    TYPE(force_env_type), POINTER            :: force_env
    LOGICAL, INTENT(IN), OPTIONAL            :: calc_force, &
                                                consistent_energies, &
                                                skip_external_control, &
                                                eval_energy_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_calc_energy_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nat, output_unit, &
                                                print_forces, stat
    LOGICAL                                  :: energy_consistency, failure, &
                                                ionode, my_eval, my_skip, &
                                                need_f, pvd, pvn
    REAL(kind=dp)                            :: e_pot
    REAL(kind=dp), DIMENSION(:), POINTER     :: pos
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_result_type), POINTER            :: results
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(virial_type), POINTER               :: virial

    NULLIFY(logger,results)
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    output_unit = -1
    IF (ionode) output_unit= cp_logger_get_default_unit_nr(logger)
    failure=.FALSE.
    need_f=.TRUE.
    my_eval=.TRUE.
    my_skip=.FALSE.
    energy_consistency=.FALSE.
    IF (PRESENT(eval_energy_forces)) my_eval = eval_energy_forces
    IF (PRESENT(skip_external_control)) my_skip = skip_external_control
    IF (PRESENT(calc_force)) need_f=calc_force
    IF (PRESENT(consistent_energies)) energy_consistency = consistent_energies

    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF
    CALL force_env_get(force_env,virial=virial,error=error)
    IF (virial%pv_availability) THEN
       pvn = virial%pv_numer
       pvd = virial%pv_diagonal
       CALL zero_virial(virial)
       CALL virial_set(virial,pv_availability=.TRUE.,pv_numer=pvn,pv_diagonal=pvd)
    END IF
    IF (.NOT.failure) THEN
       IF (my_eval) THEN 
         SELECT CASE ( force_env%in_use )
         CASE ( use_fist_force )
            CALL fist_force_control( force_env%fist_env, virial, force_env%para_env, &
                 force_env_section=force_env%force_env_section, error=error)
            CALL get_fist_env(fist_env=force_env%fist_env,results=results,error=error)
         CASE (use_ep_force)
            CALL ep_env_calc_e_f(force_env%ep_env,need_f,error=error)
         CASE ( use_qs_force )
            IF (.NOT.need_f) THEN
               CALL qs_energies(qs_env=force_env%qs_env, globenv=force_env%globenv, &
                    consistent_energies=energy_consistency, calc_forces=calc_force, error=error)
            ELSE
               CALL qs_forces(force_env%qs_env,force_env%globenv, error=error)
            END IF
            IF (virial%pv_availability.AND.virial%pv_numer) THEN
               CALL qs_calc_numerical_pressure(force_env,error=error)
            END IF
            CALL get_qs_env(qs_env=force_env%qs_env,results=results,error=error)         
         CASE ( use_kg_force )
            IF (.NOT.need_f) THEN
               root_section=>force_env%root_section
               CALL kg_energies ( force_env%kg_env,root_section, force_env%globenv, error=error )
            ELSE
               CALL kg_force_control ( force_env%kg_env, force_env%globenv, error=error )
            END IF
         CASE (use_eip_force)
            IF (force_env%eip_env%eip_model == use_lenosky_eip) THEN
               CALL eip_lenosky(force_env, error=error)
            ELSE IF (force_env%eip_env%eip_model == use_bazant_eip) THEN
               CALL eip_bazant(force_env, error=error)
            END IF
         CASE ( use_qmmm )
            CALL qmmm_energy_and_forces(force_env,need_f,error=error)
         CASE ( use_mixed_force )
            CALL mixed_energy_forces(force_env,need_f,error=error)
         CASE default
            CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
         END SELECT
       END IF
       IF (.NOT.my_skip) THEN
          ! Flexible Partitioning
          IF (ASSOCIATED(force_env%fp_env)) THEN
             IF (force_env%fp_env%use_fp) THEN
                CALL force_env_get(force_env,cell=cell,error=error)
                CALL fp_eval(force_env%fp_env,force_env%subsys, &
                     cell,error=error)
             ENDIF
          ENDIF
          ! Constraints ONLY of Fixed Atom type
          CALL fix_atom_control(force_env, error=error)
          ! All Restraints
          CALL restraint_control(force_env, error=error)
          ! Rescale forces if requested
          CALL rescale_forces(force_env, error=error)
       END IF
       IF (ASSOCIATED(force_env%cpot_env)) THEN
          CALL force_env_get(force_env, potential_energy=e_pot, error=error)
          nat=force_env_get_natom(force_env,error=error)
          ALLOCATE(pos(3*nat),stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          CALL force_env_pos_get(force_env, pos, 3*nat, error=error)
          DEALLOCATE(pos,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          CALL cpot_calc(force_env%cpot_env, e_pot,pos,error)
       END IF

       print_forces = cp_print_key_unit_nr(logger,force_env%force_env_section,"PRINT%FORCES",&
                        extension=".xyz",error=error)
       IF (print_forces>0 .AND. need_f) THEN
          CALL force_env_get(force_env, potential_energy=e_pot, error=error)
          WRITE(output_unit,'(/,T2,"ENERGY| Total FORCE_EVAL ( ",A," ) energy (a.u.): ",T55,F26.15,/)')&
               ADJUSTR(TRIM(use_prog_name(force_env%in_use))),e_pot
          CALL force_env_get(force_env, subsys=subsys, error=error)         
          CALL cp_subsys_get(subsys, particles=particles, error=error)
          CALL write_particle_coordinates(particles%els,print_forces,&
               dump_xmol,"FORCE"," FORCES| (a.u.) " )
       ENDIF
       CALL cp_print_key_finished_output(print_forces,logger,force_env%force_env_section,&
            "PRINT%FORCES",error=error)

       ! Dump force_eval energy..
       IF (output_unit>0.AND.(.NOT.need_f)) THEN
          CALL force_env_get(force_env, potential_energy=e_pot, error=error)
          WRITE(output_unit,'(/,T2,"ENERGY| Total FORCE_EVAL ( ",A," ) energy (a.u.): ",T55,F26.15,/)')&
               ADJUSTR(TRIM(use_prog_name(force_env%in_use))),e_pot
       END IF
       IF(ASSOCIATED(results))THEN
          CALL erase_result(results=force_env%results,error=error)
          CALL transfer_result_type(res_target=force_env%results,res_source=results,error=error)
       END IF
    END IF
  END SUBROUTINE force_env_calc_energy_force
!******************************************************************************

!!****f* force_env_types/force_env_create [1.0] *
!!
!!   NAME
!!     force_env_create
!!
!!   FUNCTION
!!     creates and initializes a force environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - force_env: the force env to create
!!     - fist_env, kg_env, qs_env: exactly one of these should be
!!       associated, the one that is active
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE force_env_create(force_env,root_section,para_env,globenv,fist_env,qs_env,kg_env,&
     meta_env,sub_force_env,qmmm_env,eip_env,ep_env,force_env_section,mixed_env,error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(fist_environment_type), OPTIONAL, &
      POINTER                                :: fist_env
    TYPE(qs_environment_type), OPTIONAL, &
      POINTER                                :: qs_env
    TYPE(kg_environment_type), OPTIONAL, &
      POINTER                                :: kg_env
    TYPE(meta_env_type), OPTIONAL, POINTER   :: meta_env
    TYPE(force_env_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: sub_force_env
    TYPE(qmmm_env_qm_type), OPTIONAL, &
      POINTER                                :: qmmm_env
    TYPE(eip_environment_type), OPTIONAL, &
      POINTER                                :: eip_env
    TYPE(ep_env_type), OPTIONAL, POINTER     :: ep_env
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(mixed_environment_type), OPTIONAL, &
      POINTER                                :: mixed_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat, stress_tensor
    LOGICAL                                  :: failure, pv_availability, &
                                                pv_diagonal, pv_numerical, &
                                                virial_in_use
    TYPE(cp_subsys_type), POINTER            :: subsys

  failure=.FALSE.

  ALLOCATE ( force_env, stat=stat )
  CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
  IF (.NOT. failure) THEN
     NULLIFY ( force_env%subsys, force_env%fist_env, &
          force_env%kg_env, force_env%qs_env,   &
          force_env%para_env, force_env%globenv, &
          force_env%dyn_coeff_set,&
          force_env%meta_env, force_env%sub_force_env, &
          force_env%qmmm_env, force_env%ep_env, force_env%fp_env, &
          force_env%force_env_section, force_env%eip_env,force_env%mixed_env,&
          force_env%extended_energies,force_env%root_section,&
          force_env%cpot_env,force_env%results)
     last_force_env_id=last_force_env_id+1
     force_env%id_nr=last_force_env_id
     force_env%ref_count=1
     force_env%in_use=0

     force_env%globenv => globenv
     CALL globenv_retain(force_env%globenv,error=error)

     force_env%root_section => root_section
     CALL section_vals_retain(root_section,error=error)

     force_env%para_env=>para_env
     CALL cp_para_env_retain(force_env%para_env, error=error)

     CALL section_vals_retain(force_env_section,error=error)
     force_env%force_env_section => force_env_section

     ! Should we compute the virial?
     CALL section_vals_val_get(force_env_section,"STRESS_TENSOR",i_val=stress_tensor,error=error)
     SELECT CASE(stress_tensor)
     CASE(do_stress_none)
        pv_availability=.FALSE.        
        pv_numerical=.FALSE.
        pv_diagonal=.FALSE.
     CASE(do_stress_analytical)
        pv_availability=.TRUE.        
        pv_numerical=.FALSE.
        pv_diagonal=.FALSE.
     CASE(do_stress_numerical)
        pv_availability=.TRUE.        
        pv_numerical=.TRUE.
        pv_diagonal=.FALSE.
     CASE(do_stress_diagonal_anal)
        pv_availability=.TRUE.        
        pv_numerical=.FALSE.
        pv_diagonal=.TRUE.
     CASE(do_stress_diagonal_numer)
        pv_availability=.TRUE.        
        pv_numerical=.TRUE.
        pv_diagonal=.TRUE.
     END SELECT

     IF (PRESENT(fist_env)) THEN
        IF (ASSOCIATED(fist_env)) THEN
           CPPrecondition(force_env%in_use==0,cp_failure_level,routineP,error,failure)
           force_env%in_use=use_fist_force
           force_env%fist_env => fist_env
           CALL fist_env_retain ( fist_env, error = error )
           ! Virial controlled through the external request
           CALL virial_create ( force_env % virial, error=error)
           CALL virial_set ( force_env % virial, pv_availability=pv_availability)
        END IF
     END IF
     IF (PRESENT(kg_env)) THEN
        IF (ASSOCIATED(kg_env)) THEN
           CPPrecondition(force_env%in_use==0,cp_failure_level,routineP,error,failure)
           force_env%in_use=use_kg_force
           force_env%kg_env => kg_env
           CALL kg_env_retain(kg_env,error=error)
           ! Virial not present for KG
           CALL virial_create ( force_env%virial, error=error)
           kg_env%virial => force_env%virial
           CALL virial_retain(kg_env%virial,error=error)
        END IF
     END IF
     IF (PRESENT(eip_env)) THEN
        IF (ASSOCIATED(eip_env)) THEN
           CPPrecondition(force_env%in_use==0, cp_failure_level, routineP, error, failure)
           force_env%in_use = use_eip_force
           force_env%eip_env => eip_env
           CALL eip_env_retain(eip_env, error=error)
           ! Virial not present for EIP
           CALL virial_create(force_env%virial, error=error)
           eip_env%virial => force_env%virial
           CALL virial_retain(eip_env%virial,error=error)
        END IF
     END IF
     IF (PRESENT(qs_env)) THEN
        IF (ASSOCIATED(qs_env)) THEN
           CPPrecondition(force_env%in_use==0,cp_failure_level,routineP,error,failure)
           force_env%in_use=use_qs_force
           force_env%qs_env => qs_env
           CALL qs_env_retain(qs_env,error=error)
           IF (force_env%globenv%run_type_id == qs_debug_run) THEN
              CALL section_vals_val_get(root_section,"DEBUG%DEBUG_PRESSURE",&
                   l_val=virial_in_use,error=error)
              CALL virial_create(force_env%virial, error=error)
              CALL virial_set ( force_env % virial, pv_availability = virial_in_use)
           ELSE
              CALL virial_create(force_env%virial, error=error)
              ! Virial controlled through the external request
              CALL virial_set ( virial=force_env%virial,&
                                pv_availability=pv_availability,&
                                pv_numer=pv_numerical,&
                                pv_diagonal=pv_diagonal)
           END IF
           qs_env%virial => force_env%virial
           CALL virial_retain(qs_env%virial,error=error)
        END IF
     END IF
     IF (PRESENT(qmmm_env)) THEN
        CPPrecondition(PRESENT(sub_force_env),cp_failure_level,routineP,error,failure)
        force_env%in_use=use_qmmm
        force_env%qmmm_env => qmmm_env
        CALL qmmm_env_qm_retain(qmmm_env,error=error)
        force_env%virial => sub_force_env(1)%force_env%virial
        CALL virial_retain(force_env%virial,error=error)
        ! Virial controlled through the external request
        CALL virial_set ( virial=force_env%virial,&
                          pv_availability=pv_availability,&
                          pv_numer=pv_numerical,&
                          pv_diagonal=pv_diagonal)        
     END IF
     IF (PRESENT(mixed_env)) THEN
        CPPrecondition(force_env%in_use==0, cp_failure_level, routineP, error, failure)
        force_env%in_use=use_mixed_force
        force_env%mixed_env => mixed_env
        CALL mixed_env_retain ( mixed_env, error = error )
        ! This is necessary as long as there are methods not implementing the virial
        CALL virial_create ( force_env % virial, error=error)
        CALL virial_set ( virial=force_env%virial,&
                          pv_availability=pv_availability,&
                          pv_numer=pv_numerical,&
                          pv_diagonal=pv_diagonal)
     END IF
     IF (PRESENT(ep_env)) THEN
        IF (ASSOCIATED(ep_env)) THEN
           CPPrecondition(force_env%in_use==0,cp_failure_level,routineP,error,failure)
           force_env%in_use=use_ep_force
           force_env%ep_env => ep_env
           CALL ep_env_retain(ep_env,error=error)
           ! Virial not present for EP
           CALL virial_create ( force_env%virial, error=error)
        END IF
     END IF
     CPPostcondition(force_env%in_use/=0,cp_failure_level,routineP,error,failure)

     IF (PRESENT(sub_force_env)) THEN
        force_env%sub_force_env => sub_force_env
     END IF

     IF (PRESENT(meta_env)) THEN
        force_env%meta_env => meta_env
        CALL meta_env_retain(meta_env,error=error)
     END IF

     CALL cp_result_create(results=force_env%results,error=error)
     SELECT CASE(force_env%in_use)
     CASE(use_fist_force)
        CALL get_fist_env (force_env%fist_env, subsys=force_env%subsys, error=error)
        CALL cp_subsys_retain (force_env%subsys, error=error)
     CASE(use_kg_force)
        CALL get_kg_env(force_env%kg_env, subsys=force_env%subsys,error=error)
        CALL cp_subsys_retain(force_env%subsys,error=error)
        IF(kg_env%dft_control%qs_control%polarization) THEN
           CALL get_kg_env(force_env%kg_env, dyn_coeff_set=force_env%dyn_coeff_set,&
                error=error)
           CALL dyn_coeff_set_retain(force_env%dyn_coeff_set,error=error)
        END IF
     CASE(use_qs_force)
        CALL get_qs_env(force_env%qs_env, subsys=force_env%subsys,error=error)
        CALL cp_subsys_retain(force_env%subsys,error=error)
     CASE(use_ep_force)
        CALL get_qs_env(force_env%ep_env%main_qs_env,subsys=force_env%subsys,error=error)
        CALL cp_subsys_retain(force_env%subsys,error=error)
     CASE(use_eip_force)
        CALL eip_env_get(force_env%eip_env,subsys=force_env%subsys,error=error)
        CALL cp_subsys_retain(force_env%subsys, error=error)
     CASE(use_qmmm)
        subsys => force_env%sub_force_env(1)%force_env%subsys
        force_env%subsys => subsys
        CALL cp_subsys_retain(subsys,error=error)
     CASE(use_mixed_force)
        CALL get_mixed_env (force_env%mixed_env,subsys=force_env%subsys,error=error)
        CALL cp_subsys_retain (force_env%subsys, error=error)
     CASE default
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
     END SELECT

  END IF

END SUBROUTINE force_env_create

!***************************************************************************


!!****f* qs_ep_methods/ep_create_force_env *
!!
!!   NAME
!!     ep_create_force_env
!!
!!   FUNCTION
!!     creates a force environment that does an ep calculation
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - force_env: the force environment to be created
!!     - globenv: the global environment with input,...
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE ep_create_force_env(force_env, root_section, para_env, globenv,&
     force_env_section, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_create_force_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(ep_env_type), POINTER               :: ep_env

  failure=.FALSE.

  IF (.NOT. failure) THEN
     NULLIFY(ep_env)
     CALL ep_env_create(ep_env, root_section, para_env, globenv=globenv,&
          error=error)
     CALL force_env_create(force_env,root_section,para_env,globenv=globenv,ep_env=ep_env,&
          force_env_section = force_env_section, error=error)
     CALL ep_env_release(ep_env,error=error)
  END IF
END SUBROUTINE ep_create_force_env
!***************************************************************************

!!****f* force_env_types/force_env_set_cell [1.0] *
!!
!!   NAME
!!     force_env_set_cell
!!
!!   FUNCTION
!!     changes the cell of the system
!!
!!   NOTES
!!     for future use, each method should make the changes that it needs
!!     to cope with the new cell
!!
!!   INPUTS
!!     - force_env: the force environment where the cell should be changed
!!     - cell: the new cell
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2003 created [fawzi]
!!
!!*** **********************************************************************
RECURSIVE SUBROUTINE force_env_set_cell(force_env, cell, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_set_cell', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iforce_eval, nforce_eval
    LOGICAL                                  :: failure
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(kg_rho_type), POINTER               :: kg_rho
    TYPE(kg_rspw_type), POINTER              :: kg_rspw
    TYPE(qs_ks_env_type), POINTER            :: new_ks_env

  failure=.FALSE.
  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(cell),cp_failure_level,routineP,error,failure)
  CPPrecondition(cell%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     SELECT CASE(force_env%in_use)
     CASE (use_qs_force)
        NULLIFY ( new_ks_env )
        CPPrecondition(ASSOCIATED(force_env%qs_env),cp_failure_level,routineP,error,failure)
        CALL set_qs_env(force_env%qs_env,&
             cell=cell,&
             error=error)
        CALL qs_env_rebuild_pw_env( force_env % qs_env, error=error)
        CALL qs_ks_create ( new_ks_env,  force_env % qs_env ,error=error)
        CALL set_qs_env( force_env % qs_env, ks_env=new_ks_env ,error=error)
        CALL qs_ks_release ( new_ks_env ,error=error)
     CASE (use_ep_force)
        NULLIFY ( new_ks_env )
        CPPrecondition(ASSOCIATED(force_env%ep_env),cp_failure_level,routineP,error,failure)
        CALL set_qs_env(force_env%ep_env%main_qs_env,&
             cell=cell,&
             error=error)
        CALL qs_env_rebuild_pw_env( force_env % ep_env%main_qs_env, error=error)
        CALL qs_ks_create ( new_ks_env,  force_env % qs_env ,error=error)
        CALL set_qs_env( force_env % ep_env%main_qs_env, ks_env=new_ks_env ,error=error)
        CALL qs_ks_release ( new_ks_env ,error=error)
     CASE (use_fist_force)
        CPPrecondition(ASSOCIATED(force_env%fist_env),cp_failure_level,routineP,error,failure)
        CALL get_fist_env ( force_env%fist_env, ewald_pw = ewald_pw,  &
                            ewald_env = ewald_env ,error=error)
        CALL set_fist_env(force_env%fist_env, cell=cell, error=error)
        CALL ewald_pw_grid_change ( ewald_pw, ewald_env, cell, error )
     CASE (use_kg_force)
        CPPrecondition(ASSOCIATED(force_env%kg_env),cp_failure_level,routineP,error,failure)
        CALL set_kg_env(force_env%kg_env,cell=cell, error=error)
        CALL get_kg_env ( force_env%kg_env, rspw = kg_rspw, rho=kg_rho ,error=error)
        CALL kg_rspw_rebuild( kg_rspw, force_env % kg_env, error=error)
        ! allocate rho structures
        CALL kg_rho_create(kg_rho=kg_rho,kg_env=force_env % kg_env,error=error)
        CALL set_kg_env(force_env % kg_env,rho=kg_rho,error=error)
     CASE (use_eip_force)
        CPAssert(.FALSE., cp_failure_level, routineP, error, failure)
     CASE (use_mixed_force)
        CPPrecondition(ASSOCIATED(force_env%mixed_env),cp_failure_level,routineP,error,failure)
        nforce_eval = SIZE(force_env%sub_force_env)
        DO iforce_eval = 1, nforce_eval
           IF (ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) THEN
              CALL force_env_set_cell(force_env=force_env%sub_force_env(iforce_eval)%force_env,&
                   cell=cell, error=error) 
           END IF
        END DO
     CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
             routineP,"unknown in_use flag value "//&
CPSourceFileRef,&
             error,failure)
     END SELECT
  END IF
END SUBROUTINE force_env_set_cell
!***************************************************************************

!!****f* force_env_methods/qmmm_energy_and_forces [1.0] *
!!
!!   NAME
!!     qmmm_energy_and_forces
!!
!!   FUNCTION
!!     calculates the qm/mm energy and forces
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - calc_force: if also the forces should be calculated
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [fawzi]
!!
!!*** **********************************************************************
RECURSIVE SUBROUTINE qmmm_energy_and_forces(force_env,calc_force,error)
    TYPE(force_env_type), POINTER            :: force_env
    LOGICAL, INTENT(IN), OPTIONAL            :: calc_force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_energy_and_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip, output_unit
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index
    LOGICAL                                  :: failure, need_f, &
                                                qmmm_added_chrg, qmmm_link, &
                                                qmmm_link_imomm
    REAL(KIND=dp)                            :: energy_mm, energy_qm
    REAL(KIND=dp), DIMENSION(3)              :: max_coord, min_coord
    TYPE(cell_type), POINTER                 :: mm_cell, qm_cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys_mm, subsys_qm
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_mm, particles_qm
    TYPE(qmmm_links_type), POINTER           :: qmmm_links
    TYPE(section_vals_type), POINTER         :: force_env_section

  min_coord       =  HUGE(0.0_dp)
  max_coord       = -HUGE(0.0_dp)
  failure         = .FALSE.
  need_f          = .TRUE.
  qmmm_link       = .FALSE.
  qmmm_link_imomm = .FALSE.
  qmmm_added_chrg = .FALSE.
  logger => cp_error_get_logger(error)
  IF (PRESENT(calc_force)) need_f=calc_force
  NULLIFY(subsys_mm, subsys_qm, qm_atom_index,particles_mm,particles_qm, qm_cell, mm_cell)
  NULLIFY(force_env_section)
  force_env_section => force_env%sub_force_env(qs_subsys)%force_env%force_env_section
  
  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(force_env%qmmm_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(force_env%qmmm_env%ref_count>0,cp_failure_level,routineP,error,failure)

  CALL force_env_get(force_env%sub_force_env(fist_subsys)%force_env,&
                     cell=mm_cell,subsys=subsys_mm,error=error)
  CALL force_env_get(force_env%sub_force_env(qs_subsys)%force_env,&
                     cell=qm_cell,subsys=subsys_qm,error=error)
  qm_atom_index   => force_env%qmmm_env%qm_atom_index
  qmmm_link       =  force_env%qmmm_env%qmmm_link
  qmmm_links      => force_env%qmmm_env%qmmm_links
  qmmm_added_chrg =  (force_env%qmmm_env%move_mm_charges .OR. force_env%qmmm_env%add_mm_charges)
  IF (qmmm_link) THEN
     CPPrecondition(ASSOCIATED(qmmm_links),cp_failure_level,routineP,error,failure)
     IF (ASSOCIATED(qmmm_links%imomm)) qmmm_link_imomm = (SIZE(qmmm_links%imomm) /= 0)
  END IF
  CPPrecondition(ASSOCIATED(qm_atom_index),cp_failure_level,routineP,error,failure)

  ! Possibly translate the system
  CALL apply_qmmm_translate(force_env, error)

  particles_mm => subsys_mm%particles%els
  particles_qm => subsys_qm%particles%els
  ! If present QM/MM links (just IMOMM) correct the position of the qm-link atom
  IF (qmmm_link_imomm) CALL qmmm_link_Imomm_coord(qmmm_links, particles_qm, qm_atom_index, error)
  ! If add charges get their position NOW!
  IF (qmmm_added_chrg) CALL qmmm_added_chrg_coord(force_env%qmmm_env, particles_mm, error)

  ! Initialize ks_qmmm_env
  CALL ks_qmmm_env_rebuild(qs_env=force_env%sub_force_env(qs_subsys)%force_env%qs_env,&
       qmmm_env=force_env%qmmm_env,error=error)

  ! Compute the short range QM/MM Electrostatic Potential
  CALL qmmm_el_coupling( qs_env=force_env%sub_force_env(qs_subsys)%force_env%qs_env,&
       qmmm_env=force_env%qmmm_env,&
       mm_particles=particles_mm,&
       mm_cell=mm_cell,&
       error=error)

  ! Fist
  CALL force_env_calc_energy_force(force_env%sub_force_env(fist_subsys)%force_env,&
       calc_force=need_f,skip_external_control=.TRUE.,error=error)
  !
  ! Print Out information on fist energy calculation...
  !
  CALL force_env_get(force_env%sub_force_env(fist_subsys)%force_env,&
                     potential_energy=energy_mm,&
                     error=error)
  !
  ! QS
  CALL force_env_calc_energy_force(force_env%sub_force_env(qs_subsys)%force_env,&
       calc_force=need_f,skip_external_control=.TRUE.,error=error)
  !
  ! Print Out information on QS energy calculation...
  !
  CALL force_env_get(force_env%sub_force_env(qs_subsys)%force_env,&
                     potential_energy=energy_qm,&
                     error=error)
  !
  ! QM/MM Interaction Potential forces
  !
  CALL qmmm_forces(force_env%sub_force_env(qs_subsys)%force_env%qs_env,&
       force_env%qmmm_env,particles_mm,&
       mm_cell=mm_cell,&
       calc_force=need_f,error=error)
  !
  ! Print Out information on QS energy calculation...
  !
  CALL force_env_get(force_env%sub_force_env(qs_subsys)%force_env,&
                     potential_energy=energy_qm,&
                     error=error)

  IF (need_f) THEN
     ! If present QM/MM links (just IMOMM) correct the position of the qm-link atom
     IF (qmmm_link_imomm) CALL qmmm_link_Imomm_forces(qmmm_links,particles_qm,qm_atom_index,error)
     particles_mm => subsys_mm%particles%els
     DO ip=1,SIZE(qm_atom_index)
        particles_mm(qm_atom_index(ip))%f=particles_mm(qm_atom_index(ip))%f+particles_qm(ip)%f
     END DO
     ! If add charges get rid of their derivatives right NOW!
     IF (qmmm_added_chrg) CALL qmmm_added_chrg_forces(force_env%qmmm_env, particles_mm, error)
  END IF

  output_unit = cp_print_key_unit_nr(logger,force_env_section,"QMMM%PRINT%DERIVATIVES",&
           extension=".Log",error=error)
  IF (output_unit>0) THEN
     WRITE (unit=output_unit,fmt='(/1X,A,F15.9)')"Energy after QMMM calculation: ",energy_qm
     IF (need_f) THEN
        WRITE (unit=output_unit,fmt='(/1X,A)')"Derivatives on all atoms after QMMM calculation: "
        DO ip=1,SIZE(particles_mm)
           WRITE (unit=output_unit,fmt='(1X,3F15.9)')particles_mm(ip)%f
        END DO
     END IF
  END IF
  CALL cp_print_key_finished_output(output_unit,logger,force_env_section,&
       "QMMM%PRINT%DERIVATIVES",error=error)
END SUBROUTINE qmmm_energy_and_forces

!! ****f* force_env_methods/mixed_energy_forces  [1.0] 
!!
!!   NAME
!!    mixed_energy_forces
!!
!!   FUNCTION
!!    Computes energy and forces for a mixed force_env type
!!
!!   NOTES
!!    -
!!
!!   INPUTS
!!
!!
!!   AUTHOR
!!     Florian Schiffmann
!!
!!   MODIFICATION HISTORY
!!      11.06  created [fschiff]
!!      04.07  generalization to an illimited number of force_eval [tlaino]
!!      04.07  further generalization to force_eval with different geometrical
!!             structures
!!
!!***************************************************************************
SUBROUTINE mixed_energy_forces(force_env, need_f, error)
    TYPE(force_env_type), POINTER            :: force_env
    LOGICAL, INTENT(IN)                      :: need_f
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mixed_energy_forces', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: coupling_function
    CHARACTER(LEN=default_string_length)     :: def_error, this_error
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: data_function, my_par, &
                                                parameters, variables
    INTEGER :: iforce_eval, iparticle, j, jparticle, mixing_type, my_group, &
      mydim, natom, nforce_eval, output_unit, stat
    INTEGER, DIMENSION(:), POINTER           :: glob_natoms, map_index
    LOGICAL                                  :: failure, ionode, &
                                                virial_consistent
    REAL(KIND=dp) :: coupling_parameter, dedf, der_1, der_2, dx, energy, err, &
      lambda, lerr, restraint_strength, restraint_target, sd
    REAL(KIND=dp), DIMENSION(:), POINTER     :: energies, my_val, par_values
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: forces
    TYPE(cell_type), POINTER                 :: cell, cell_mix
    TYPE(cp_error_type)                      :: my_error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_p_type), DIMENSION(:), &
      POINTER                                :: subsystems
    TYPE(cp_subsys_type), POINTER            :: subsys_mix
    TYPE(mixed_energy_type), POINTER         :: mixed_energy
    TYPE(mixed_force_type), DIMENSION(:), &
      POINTER                                :: global_forces
    TYPE(particle_list_p_type), &
      DIMENSION(:), POINTER                  :: particles
    TYPE(particle_list_type), POINTER        :: particles_mix
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                mapping_section, root_section
    TYPE(virial_p_type), DIMENSION(:), &
      POINTER                                :: virials
    TYPE(virial_type), POINTER               :: loc_virial, virial_mix

    failure=.FALSE.
    output_unit = -1
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) output_unit= cp_logger_get_default_unit_nr(logger)    
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
    ! Get infos about the mixed subsys
    CALL force_env_get(force_env=force_env,&
                       subsys=subsys_mix,&
                       force_env_section=force_env_section,&
                       root_section=root_section,&
                       virial=virial_mix,&
                       cell=cell_mix,&
                       error=error)
    CALL cp_subsys_get(subsys=subsys_mix,&
                       particles=particles_mix,&
                       error=error)
    NULLIFY(map_index, glob_natoms, global_forces)
    virial_consistent = .TRUE.
    nforce_eval = SIZE(force_env%sub_force_env)
    mapping_section => section_vals_get_subs_vals(force_env_section,"MIXED%MAPPING",error=error)
    ! Global Info
    ALLOCATE(subsystems(nforce_eval), stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(particles(nforce_eval), stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ! Local Info to sync
    ALLOCATE(global_forces(nforce_eval), stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(energies(nforce_eval), stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(glob_natoms(nforce_eval), stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(virials(nforce_eval), stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    energies    = 0.0_dp
    glob_natoms = 0
    DO iforce_eval = 1, nforce_eval
       NULLIFY(subsystems(iforce_eval)%subsys, particles(iforce_eval)%list, virials(iforce_eval)%virial)
       CALL virial_create (virials(iforce_eval)%virial, error)
       IF (.NOT.ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) CYCLE
       ! From this point on the error is the sub_error
       my_group = force_env%mixed_env%group_distribution(force_env%para_env%mepos)
       my_error = force_env%mixed_env%sub_error(my_group+1)
       ! Get all available subsys
       CALL force_env_get(force_env=force_env%sub_force_env(iforce_eval)%force_env,&
                          subsys=subsystems(iforce_eval)%subsys,cell=cell,error=my_error)
       ! Check whether virial can be consistently used..
       IF (virial_mix%pv_availability) THEN
          virial_consistent = virial_consistent.AND.compare_cells(cell_mix, cell, my_error)
       END IF
       ! Get available particles
       CALL cp_subsys_get(subsys=subsystems(iforce_eval)%subsys,&
                          particles=particles(iforce_eval)%list,error=my_error)

       ! Get Mapping index array
       natom = SIZE(particles(iforce_eval)%list%els)
       CALL get_subsys_map_index(mapping_section, natom, iforce_eval, nforce_eval, &
            map_index, my_error)

       ! Mapping particles from iforce_eval environment to the mixed env
       DO iparticle = 1, natom
          jparticle = map_index(iparticle)
          particles(iforce_eval)%list%els(iparticle)%r= particles_mix%els(jparticle)%r
       END DO

       ! Calculate energy and forces for each sub_force_env
       CALL force_env_calc_energy_force(force_env%sub_force_env(iforce_eval)%force_env,&
                                        calc_force=need_f,&
                                        skip_external_control=.TRUE.,&
                                        error=my_error)
       ! Only the rank 0 process collect info for each computation
       IF ( force_env%sub_force_env(iforce_eval)%force_env%para_env%mepos==&
            force_env%sub_force_env(iforce_eval)%force_env%para_env%source) THEN
          CALL force_env_get(force_env%sub_force_env(iforce_eval)%force_env,&
                             potential_energy=energy,&
                             virial=loc_virial,&
                             error=my_error)
          energies(iforce_eval)    = energy
          glob_natoms(iforce_eval) = natom 
          CALL cp_virial (loc_virial, virials(iforce_eval)%virial)
       END IF
       ! Deallocate map_index array
       IF (ASSOCIATED(map_index)) THEN
          DEALLOCATE(map_index, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,my_error,failure)
       END IF
       CALL cp_error_check(my_error, failure)
    END DO
    ! Final check on virial
    CALL cp_assert(virial_consistent,cp_failure_level,cp_assertion_failed,&
         routineP,"Mixed force_eval have different cells definition. Virial cannot be "//&
         " defined in a consistent way. Check the CELL sections! "//&
 CPSourceFileRef,&
         error,failure)

    ! Handling Parallel execution 
    CALL mp_sync(force_env%para_env%group)
    ! Let's transfer energy, natom, forces, virials
    CALL mp_sum(energies, force_env%para_env%group)
    CALL mp_sum(glob_natoms, force_env%para_env%group)
    ! Transfer forces
    DO iforce_eval = 1, nforce_eval
       ALLOCATE(global_forces(iforce_eval)%forces(3,glob_natoms(iforce_eval)),stat=stat)
       global_forces(iforce_eval)%forces = 0.0_dp
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) THEN
          IF ( force_env%sub_force_env(iforce_eval)%force_env%para_env%mepos==&
               force_env%sub_force_env(iforce_eval)%force_env%para_env%source) THEN
             ! Forces
             DO iparticle = 1, glob_natoms(iforce_eval)
                global_forces(iforce_eval)%forces(:,iparticle) = &
                     particles(iforce_eval)%list%els(iparticle)%f
             END DO
          END IF
       END IF
       CALL mp_sum(global_forces(iforce_eval)%forces, force_env%para_env%group)
       !Transfer only the relevant part of the virial..
       CALL mp_sum(virials(iforce_eval)%virial%pv_total, force_env%para_env%group)
       CALL mp_sum(virials(iforce_eval)%virial%pv_kinetic, force_env%para_env%group)
       CALL mp_sum(virials(iforce_eval)%virial%pv_virial, force_env%para_env%group)
       CALL mp_sum(virials(iforce_eval)%virial%pv_xc, force_env%para_env%group)
       CALL mp_sum(virials(iforce_eval)%virial%pv_constraint, force_env%para_env%group)
    END DO
    force_env%mixed_env%energies = energies
    ! Start combining the different sub_force_env
    CALL get_mixed_env(mixed_env=force_env%mixed_env,&
                       mixed_energy=mixed_energy,&
                       error=error)

    CALL section_vals_val_get(force_env_section,"MIXED%MIXING_TYPE",i_val=mixing_type,error=error)
    SELECT CASE(mixing_type)
    CASE(mix_linear_combination)
       ! Support offered only 2 force_eval
       CPPrecondition(nforce_eval==2,cp_failure_level,routineP,error,failure)
       CALL section_vals_val_get(force_env_section,"MIXED%LINEAR%LAMBDA",&
            r_val=lambda,error=error)
       mixed_energy%pot=lambda*energies(1) + (1.0_dp-lambda)*energies(2)
       ! General Mapping of forces...
       CALL mixed_map_forces(particles_mix, virial_mix, global_forces, virials, lambda, 1, nforce_eval,&
            map_index, mapping_section, .TRUE., error)
       CALL mixed_map_forces(particles_mix, virial_mix, global_forces, virials, (1.0_dp-lambda), 2, nforce_eval,&
            map_index, mapping_section, .FALSE., error)
    CASE(mix_minimum)
       ! Support offered only 2 force_eval
       CPPrecondition(nforce_eval==2,cp_failure_level,routineP,error,failure)
       IF (energies(1)<energies(2)) THEN
          mixed_energy%pot=energies(1)
          CALL mixed_map_forces(particles_mix, virial_mix, global_forces, virials, 1.0_dp, 1, nforce_eval,&
               map_index, mapping_section, .TRUE., error)
       ELSE
          mixed_energy%pot=energies(2)
          CALL mixed_map_forces(particles_mix, virial_mix, global_forces, virials, 1.0_dp, 2, nforce_eval,&
               map_index, mapping_section, .TRUE., error)
       ENDIF
    CASE(mix_coupled)
       ! Support offered only 2 force_eval
       CPPrecondition(nforce_eval==2,cp_failure_level,routineP,error,failure)
       CALL section_vals_val_get(force_env_section,"MIXED%COUPLING%COUPLING_PARAMETER",&
            r_val=coupling_parameter,error=error)
       sd = SQRT((energies(1)-energies(2))**2+4.0_dp*coupling_parameter**2)
       der_1=(1.0_dp-(1.0_dp/(2.0_dp*sd))*2.0_dp*(energies(1)-energies(2)))/2.0_dp
       der_2=(1.0_dp+(1.0_dp/(2.0_dp*sd))*2.0_dp*(energies(1)-energies(2)))/2.0_dp
       mixed_energy%pot=(energies(1)+energies(2)-sd)/2.0_dp
       ! General Mapping of forces...
       CALL mixed_map_forces(particles_mix, virial_mix, global_forces, virials, der_1, 1, nforce_eval,&
            map_index, mapping_section, .TRUE., error)
       CALL mixed_map_forces(particles_mix, virial_mix, global_forces, virials, der_2, 2, nforce_eval,&
            map_index, mapping_section, .FALSE., error)
    CASE(mix_restrained)
       ! Support offered only 2 force_eval
       CPPrecondition(nforce_eval==2,cp_failure_level,routineP,error,failure)
       CALL section_vals_val_get(force_env_section,"MIXED%RESTRAINT%RESTRAINT_TARGET",&
            r_val=restraint_target,error=error)
       CALL section_vals_val_get(force_env_section,"MIXED%RESTRAINT%RESTRAINT_STRENGTH",&
            r_val=restraint_strength,error=error)
       mixed_energy%pot=energies(1)+restraint_strength*(energies(1)-energies(2)-restraint_target)**2
       der_2 = -2.0_dp*restraint_strength*(energies(1)-energies(2)-restraint_target)
       der_1 = 1.0_dp - der_2
       ! General Mapping of forces...
       CALL mixed_map_forces(particles_mix, virial_mix, global_forces, virials, der_1, 1, nforce_eval,&
            map_index, mapping_section, .TRUE., error)
       CALL mixed_map_forces(particles_mix, virial_mix, global_forces, virials, der_2, 2, nforce_eval,&
            map_index, mapping_section, .FALSE., error)
    CASE(mix_generic)
       ! Support any number of force_eval sections
       NULLIFY(data_function, variables, parameters, par_values)
       CALL section_vals_val_get(force_env_section,"MIXED%GENERIC%MIXING_FUNCTION",&
            c_vals=data_function,error=error)
       ! Clean the function structure
       coupling_function = ""
       DO j = 1, SIZE(data_function)
          coupling_function(LEN_TRIM(coupling_function)+1:) = TRIM(data_function(j))
       END DO
       CALL section_vals_val_get(force_env_section,"MIXED%GENERIC%VARIABLES",&
            c_vals=variables,error=error)
       CALL section_vals_val_get(force_env_section,"MIXED%GENERIC%PARAMETERS",&
            c_vals=parameters,error=error)
       CALL section_vals_val_get(force_env_section,"MIXED%GENERIC%VALUES",&
            r_vals=par_values,error=error)
       mydim = SIZE(variables)
       IF (ASSOCIATED(parameters)) THEN
          CPPrecondition(SIZE(parameters)==SIZE(par_values),cp_failure_level,routineP,error,failure)
          mydim=mydim+SIZE(parameters)
       END IF
       ALLOCATE(my_par(mydim),stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(my_val(mydim),stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       my_par(1:SIZE(variables)) = variables
       my_val(1:SIZE(variables)) = energies
       IF (ASSOCIATED(parameters)) THEN
          my_par(SIZE(variables)+1:) = parameters(1:SIZE(parameters))
          my_val(SIZE(variables)+1:) = par_values(1:SIZE(parameters))
       END IF
       CALL initf(1)
       CALL parsef(1,TRIM(coupling_function),my_par)
       ! Now the hardest part.. map energy with corresponding force_eval
       mixed_energy%pot= evalf(1,my_val)
       CPPrecondition(EvalErrType<=0,cp_failure_level,routineP,error,failure)
       DO iparticle = 1, SIZE(particles_mix%els)
          particles_mix%els(iparticle)%f(:) = 0.0_dp
       END DO
       CALL zero_virial(virial_mix, reset=.FALSE.)
       DO iforce_eval = 1, nforce_eval
          CALL section_vals_val_get(force_env_section,"MIXED%GENERIC%DX",r_val=dx,error=error)
          CALL section_vals_val_get(force_env_section,"MIXED%GENERIC%ERROR_LIMIT",r_val=lerr,error=error)
          dedf = evalfd(1,iforce_eval,my_val,dx,err)
          IF ((ABS(err)>lerr).AND.(output_unit>0)) THEN
             WRITE(this_error,"(A,G12.6,A)")"(",err,")"
             WRITE(def_error,"(A,G12.6,A)")"(",lerr,")"
             CALL compress(this_error,.TRUE.)
             CALL compress(def_error,.TRUE.)
             CALL cp_assert(.FALSE.,cp_warning_level,-300,routineP,&
                  'ASSERTION (cond) failed at line '//cp_to_string(__LINE__)//&
                  ' Error '//TRIM(this_error)//' in computing numerical derivatives larger then'//&
                  TRIM(def_error)//' .',error)
          END IF
          ! General Mapping of forces...
          CALL mixed_map_forces(particles_mix, virial_mix, global_forces, virials, dedf, iforce_eval,&
               nforce_eval, map_index, mapping_section, .FALSE., error)
          my_val(iforce_eval) = energies(iforce_eval)
       END DO
       ! Let's store the needed information..
       force_env%mixed_env%par = my_par
       force_env%mixed_env%val = my_val
       force_env%mixed_env%dx  = dx
       force_env%mixed_env%lerr= lerr
       force_env%mixed_env%coupling_function = TRIM(coupling_function)
       CALL finalizef()
       DEALLOCATE(my_par,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(my_val,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    CASE DEFAULT
       CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT
    !Simply deallocate and loose the pointer references..
    DO iforce_eval = 1, nforce_eval
       DEALLOCATE(global_forces(iforce_eval)%forces,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL virial_release(virials(iforce_eval)%virial, error=error)
    END DO
    DEALLOCATE(global_forces, stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(subsystems, stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(particles, stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(energies, stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(glob_natoms, stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(virials, stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
END SUBROUTINE mixed_energy_forces

!!****f* force_env_types/mixed_map_forces *
!!
!!   NAME
!!     mixed_map_forces
!!
!!   FUNCTION
!!     Maps forces between the different force_eval sections/environments
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - force_env: the force env to shake
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino - University of Zurich [tlaino] - 05.2007
!!
!!*** **********************************************************************
  SUBROUTINE mixed_map_forces(particles_mix, virial_mix, global_forces, virials, factor, iforce_eval,&
     nforce_eval, map_index, mapping_section, overwrite, error)
    TYPE(particle_list_type), POINTER        :: particles_mix
    TYPE(virial_type), POINTER               :: virial_mix
    TYPE(mixed_force_type), DIMENSION(:), &
      POINTER                                :: global_forces
    TYPE(virial_p_type), DIMENSION(:), &
      POINTER                                :: virials
    REAL(KIND=dp), INTENT(IN)                :: factor
    INTEGER, INTENT(IN)                      :: iforce_eval, nforce_eval
    INTEGER, DIMENSION(:), POINTER           :: map_index
    TYPE(section_vals_type), POINTER         :: mapping_section
    LOGICAL, INTENT(IN)                      :: overwrite
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mixed_map_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iparticle, jparticle, natom, &
                                                stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    ! Get Mapping index array
    natom = SIZE(global_forces(iforce_eval)%forces,2)
    CALL get_subsys_map_index(mapping_section, natom, iforce_eval, nforce_eval, map_index, error)
    DO iparticle = 1, natom
       jparticle = map_index(iparticle)
       IF (overwrite) THEN
          particles_mix%els(jparticle)%f(:)= factor* global_forces(iforce_eval)%forces(:,iparticle)
       ELSE
          particles_mix%els(jparticle)%f(:)= particles_mix%els(jparticle)%f(:) + &
               factor* global_forces(iforce_eval)%forces(:,iparticle)
       END IF
    END DO
    ! Mixing Virial
    IF (virial_mix%pv_availability) THEN
       IF (overwrite) CALL zero_virial(virial_mix,reset=.FALSE.)
       virial_mix%pv_total      = virial_mix%pv_total + factor*virials(iforce_eval)%virial%pv_total
       virial_mix%pv_kinetic    = virial_mix%pv_kinetic + factor*virials(iforce_eval)%virial%pv_kinetic
       virial_mix%pv_virial     = virial_mix%pv_virial + factor*virials(iforce_eval)%virial%pv_virial
       virial_mix%pv_xc         = virial_mix%pv_xc + factor*virials(iforce_eval)%virial%pv_xc
       virial_mix%pv_constraint = virial_mix%pv_constraint + factor*virials(iforce_eval)%virial%pv_constraint
    END IF
    ! Deallocate map_index array
    IF (ASSOCIATED(map_index)) THEN
       DEALLOCATE(map_index, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    
  END SUBROUTINE mixed_map_forces

!!****f* force_env_types/get_subsys_map_index *
!!
!!   NAME
!!     get_subsys_map_index
!!
!!   FUNCTION
!!     performs mapping of the subsystems of different force_eval
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - force_env: the force env to shake
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino - University of Zurich [tlaino] - 05.2007
!!
!!*** **********************************************************************
  SUBROUTINE get_subsys_map_index(mapping_section, natom, iforce_eval, nforce_eval, map_index, error)
    TYPE(section_vals_type), POINTER         :: mapping_section
    INTEGER, INTENT(IN)                      :: natom, iforce_eval, &
                                                nforce_eval
    INTEGER, DIMENSION(:), POINTER           :: map_index
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_subsys_map_index', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iatom, ival, j, jval, k, &
                                                n_rep, n_rep_loc, n_rep_map, &
                                                n_rep_sys, stat, tmp
    INTEGER, DIMENSION(:), POINTER           :: index_glo, index_loc, list
    LOGICAL                                  :: check, explicit, failure
    TYPE(section_vals_type), POINTER         :: fragments_loc, fragments_sys, &
                                                map_force_ev, map_full_sys

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(map_index),cp_failure_level,routineP,error,failure)
    ALLOCATE(map_index(natom),stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL section_vals_get(mapping_section, explicit=explicit, error=error)
    IF (.NOT.explicit) THEN
       ! Standard Mapping.. subsys are assumed to have the same structure
       DO i = 1, natom
          map_index(i) = i
       END DO
    ELSE
       ! Mapping systems with different structures
       map_full_sys => section_vals_get_subs_vals(mapping_section,"FORCE_EVAL_MIXED",error=error)
       map_force_ev => section_vals_get_subs_vals(mapping_section,"FORCE_EVAL",error=error)
       CALL section_vals_get(map_full_sys, explicit=explicit, error=error)
       CPPrecondition(explicit,cp_failure_level,routineP,error,failure)
       CALL section_vals_get(map_force_ev, explicit=explicit, n_repetition=n_rep, error=error)
       CPPrecondition(explicit,cp_failure_level,routineP,error,failure)
       CPPrecondition(n_rep==nforce_eval,cp_failure_level,routineP,error,failure)
       DO i = 1, n_rep
          CALL section_vals_val_get(map_force_ev,"_SECTION_PARAMETERS_",i_rep_section=i,i_val=ival,error=error)
          IF (ival==iforce_eval) EXIT
       END DO
       CPPrecondition(i<=nforce_eval,cp_failure_level,routineP,error,failure)
       fragments_sys => section_vals_get_subs_vals(map_full_sys,"FRAGMENT",error=error)
       fragments_loc => section_vals_get_subs_vals(map_force_ev,"FRAGMENT",i_rep_section=i,error=error)
       !Perform few check on the structure of the input mapping section. as provided by the user
       CALL section_vals_get(fragments_loc, n_repetition=n_rep_loc, error=error)
       CALL section_vals_get(fragments_sys, explicit=explicit, n_repetition=n_rep_sys, error=error)
       CPPrecondition(explicit,cp_failure_level,routineP,error,failure)
       CPPrecondition(n_rep_sys>=n_rep_loc,cp_failure_level,routineP,error,failure)
       IF (n_rep_loc==0) THEN
          NULLIFY(list)
          ! We expect an easier syntax in this case..
          CALL section_vals_val_get(map_force_ev,"DEFINE_FRAGMENTS",i_rep_section=i,n_rep_val=n_rep_map,&
               error=error)
          check = (n_rep_map/=0)
          CPPrecondition(check,cp_failure_level,routineP,error,failure)
          CALL section_vals_val_get(map_force_ev,"DEFINE_FRAGMENTS",i_rep_section=i,i_vals=list,&
               error=error)
          CPPrecondition(SIZE(list)>0,cp_failure_level,routineP,error,failure)
          iatom = 0
          DO i = 1, SIZE(list)
             jval = list(i)
             DO j = 1, n_rep_sys
                CALL section_vals_val_get(fragments_sys,"_SECTION_PARAMETERS_",i_rep_section=j,i_val=tmp,error=error)
                IF (tmp==jval) EXIT
             END DO
             CALL section_vals_val_get(fragments_sys,"_DEFAULT_KEYWORD_",i_rep_section=j,i_vals=index_glo,error=error)
             DO k = 0, index_glo(2)-index_glo(1)
                iatom = iatom + 1
 CPPrecondition(iatom<=natom,cp_failure_level,routineP,error,failure)
                map_index(iatom) = index_glo(1)+k
             END DO
          END DO
          check = (iatom==natom)
          CPPrecondition(check,cp_failure_level,routineP,error,failure)
       ELSE
          ! General syntax..
          !Loop over the fragment of the force_eval
          DO i = 1, n_rep_loc
             CALL section_vals_val_get(fragments_loc,"_SECTION_PARAMETERS_",i_rep_section=i,i_val=ival,error=error)
             CALL section_vals_val_get(fragments_loc,"MAP",i_rep_section=i,i_val=jval,error=error)
             ! Index corresponding to the mixed_force_eval fragment
             DO j = 1, n_rep_sys
                CALL section_vals_val_get(fragments_sys,"_SECTION_PARAMETERS_",i_rep_section=j,i_val=tmp,error=error)
                IF (tmp==jval) EXIT
             END DO
             CPPrecondition(j<=n_rep_sys,cp_failure_level,routineP,error,failure)
             CALL section_vals_val_get(fragments_loc,"_DEFAULT_KEYWORD_",i_rep_section=i,i_vals=index_loc,error=error)
             CALL section_vals_val_get(fragments_sys,"_DEFAULT_KEYWORD_",i_rep_section=j,i_vals=index_glo,error=error)
             check = ((index_loc(2)-index_loc(1))==(index_glo(2)-index_glo(1)))
             CPPrecondition(check,cp_failure_level,routineP,error,failure)
             ! Now let's build the real mapping
             DO k = 0, index_loc(2)-index_loc(1)
                map_index(index_loc(1)+k) = index_glo(1)+k
             END DO
          END DO
       END IF
    END IF
    
  END SUBROUTINE get_subsys_map_index
     
!!****f* force_env_types/force_env_shake *
!!
!!   NAME
!!     force_env_shake
!!
!!   FUNCTION
!!     perform shake (enforcing of constraints)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - force_env: the force env to shake
!!     - dt: the dt for shake (if you are not interested in the velocities
!!       it can be any positive number)
!!     - shake_tol: the tolerance for shake
!!     - log_unit: if >0 then some information on the shake is printed,
!!       defaults to -1
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE force_env_shake(force_env,dt,shake_tol,log_unit,lagrange_mult,dump_lm,&
       pos,vel,compold,reset,error)
    TYPE(force_env_type), POINTER            :: force_env
    REAL(kind=dp), INTENT(IN), OPTIONAL      :: dt
    REAL(kind=dp), INTENT(IN)                :: shake_tol
    INTEGER, INTENT(in), OPTIONAL            :: log_unit, lagrange_mult
    LOGICAL, INTENT(IN), OPTIONAL            :: dump_lm
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL, TARGET        :: pos, vel
    LOGICAL, INTENT(IN), OPTIONAL            :: compold, reset
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_shake', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iparticle, iparticle_kind, iparticle_local, j, &
      my_lagrange_mult, my_log_unit, nparticle_kind, nparticle_local, stat
    LOGICAL                                  :: failure, has_pos, has_vel, &
                                                my_dump_lm
    REAL(KIND=dp)                            :: mydt
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: my_pos, my_vel
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(particle_list_type), POINTER        :: particles

    failure=.FALSE.

    CALL timeset(routineN,handle)
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       my_log_unit=-1
       IF (PRESENT(log_unit)) my_log_unit=log_unit
       my_lagrange_mult=-1
       IF (PRESENT(lagrange_mult)) my_lagrange_mult=lagrange_mult
       my_dump_lm = .FALSE.
       IF (PRESENT(dump_lm)) my_dump_lm = dump_lm
       NULLIFY(subsys,cell,molecules,molecule_kinds,local_molecules,particles,&
            my_pos,my_vel,gci)
       IF (PRESENT(pos)) my_pos => pos
       IF (PRESENT(vel)) my_vel => vel
       mydt = 0.1_dp
       IF (PRESENT(dt)) mydt = dt
       CALL force_env_get(force_env,subsys=subsys,cell=cell,error=error)
       CALL cp_subsys_get(subsys, &
            atomic_kinds=atomic_kinds,&
            local_molecules_new=local_molecules,&
            local_particles=local_particles,&
            molecules_new=molecules,&
            molecule_kinds_new=molecule_kinds,&
            particles=particles,&
            gci=gci,&
            error=error)
       nparticle_kind = atomic_kinds%n_els
       IF (PRESENT(compold)) THEN
          IF (compold) THEN
             CALL getold( gci, local_molecules, molecules%els, molecule_kinds%els,&
                  particles%els, cell, error)
          END IF
       END IF
       has_pos=.FALSE.
       IF (.NOT.ASSOCIATED(my_pos)) THEN
          has_pos=.TRUE.
          ALLOCATE(my_pos(3,particles%n_els),stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          my_pos = 0.0_dp
          DO iparticle_kind=1,nparticle_kind
             nparticle_local = local_particles%n_el(iparticle_kind)
             DO iparticle_local=1,nparticle_local
                iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                my_pos (:,iparticle) = particles%els(iparticle)%r(:)
             END DO
          END DO
       END IF
       has_vel=.FALSE.
       IF (.NOT.ASSOCIATED(my_vel)) THEN
          has_vel=.TRUE.
          ALLOCATE(my_vel(3,particles%n_els),stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          my_vel = 0.0_dp
          DO iparticle_kind=1,nparticle_kind
             nparticle_local = local_particles%n_el(iparticle_kind)
             DO iparticle_local=1,nparticle_local
                iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                my_vel (:,iparticle) = particles%els(iparticle)%v(:)
             END DO
          END DO
       END IF

       CALL shake_control( gci=gci, local_molecules=local_molecules,&
            molecule_set=molecules%els, molecule_kind_set=molecule_kinds%els, &
            particle_set=particles%els, pos=my_pos, vel=my_vel, dt=mydt,&
            shake_tol=shake_tol, log_unit=my_log_unit, lagrange_mult=my_lagrange_mult,&
            dump_lm= my_dump_lm, cell=cell,group=force_env%para_env%group,&
            local_particles=local_particles, error=error )

       ! Possibly reset the lagrange multipliers
       IF (PRESENT(reset)) THEN
          IF (reset) THEN
             ! Reset Intramolecular constraints
             DO i = 1,SIZE(molecules%els)
                IF (ASSOCIATED(molecules%els(i)%lci%lcolv)) THEN
                   DO j = 1, SIZE(molecules%els(i)%lci%lcolv)
                      ! Reset langrange multiplier
                      molecules%els(i)%lci%lcolv(j)%lambda = 0.0_dp
                   END DO
                END IF
                IF (ASSOCIATED(molecules%els(i)%lci%lg3x3)) THEN
                   DO j = 1, SIZE(molecules%els(i)%lci%lg3x3)
                      ! Reset langrange multiplier
                      molecules%els(i)%lci%lg3x3(j)%lambda = 0.0_dp
                   END DO
                END IF
                IF (ASSOCIATED(molecules%els(i)%lci%lg4x6)) THEN
                   DO j = 1, SIZE(molecules%els(i)%lci%lg4x6)
                      ! Reset langrange multiplier
                      molecules%els(i)%lci%lg4x6(j)%lambda = 0.0_dp
                   END DO
                END IF
             END DO
             ! Reset Intermolecular constraints
             IF (ASSOCIATED(gci)) THEN
                IF (ASSOCIATED(gci%lcolv)) THEN
                   DO j = 1, SIZE(gci%lcolv)
                      ! Reset langrange multiplier
                      gci%lcolv(j)%lambda = 0.0_dp
                   END DO
                END IF
                IF (ASSOCIATED(gci%lg3x3)) THEN
                   DO j = 1, SIZE(gci%lg3x3)
                      ! Reset langrange multiplier
                      gci%lg3x3(j)%lambda = 0.0_dp
                   END DO
                END IF
                IF (ASSOCIATED(gci%lg4x6)) THEN
                   DO j = 1, SIZE(gci%lg4x6)
                      ! Reset langrange multiplier
                      gci%lg4x6(j)%lambda = 0.0_dp
                   END DO
                END IF
             END IF
          END IF
       END IF

       IF (has_pos) THEN
          CALL update_particle_set ( particles%els, force_env%para_env%group, pos=my_pos)
          DEALLOCATE(my_pos,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       IF (has_vel) THEN
          CALL update_particle_set ( particles%els, force_env%para_env%group, vel=my_vel)
          DEALLOCATE(my_vel,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE force_env_shake
!***************************************************************************

!!****f* force_env_types/force_env_rattle *
!!
!!   NAME
!!     force_env_rattle
!!
!!   FUNCTION
!!     perform rattle (enforcing of constraints on velocities)
!!     This routine can be easily adapted to performe rattle on whatever
!!     other vector different from forces..
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - force_env: the force env to shake
!!     - dt: the dt for shake (if you are not interested in the velocities
!!       it can be any positive number)
!!     - shake_tol: the tolerance for shake
!!     - log_unit: if >0 then some information on the shake is printed,
!!       defaults to -1
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     tlaino
!!
!!*** **********************************************************************
  SUBROUTINE force_env_rattle(force_env,dt,shake_tol,log_unit,lagrange_mult,dump_lm,&
       vel,reset,error)
    TYPE(force_env_type), POINTER            :: force_env
    REAL(kind=dp), INTENT(in), OPTIONAL      :: dt
    REAL(kind=dp), INTENT(in)                :: shake_tol
    INTEGER, INTENT(in), OPTIONAL            :: log_unit, lagrange_mult
    LOGICAL, INTENT(IN), OPTIONAL            :: dump_lm
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL, TARGET        :: vel
    LOGICAL, INTENT(IN), OPTIONAL            :: reset
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_rattle', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iparticle, iparticle_kind, iparticle_local, j, &
      my_lagrange_mult, my_log_unit, nparticle_kind, nparticle_local, stat
    LOGICAL                                  :: failure, has_vel, my_dump_lm
    REAL(KIND=dp)                            :: mydt
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: my_vel
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(particle_list_type), POINTER        :: particles

    failure=.FALSE.

    CALL timeset(routineN,handle)
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       my_log_unit=-1
       IF (PRESENT(log_unit)) my_log_unit=log_unit
       my_lagrange_mult=-1
       IF (PRESENT(lagrange_mult)) my_lagrange_mult=lagrange_mult
       my_dump_lm = .FALSE.
       IF (PRESENT(dump_lm)) my_dump_lm = dump_lm
       NULLIFY(subsys,cell,molecules,molecule_kinds,local_molecules,particles,&
            my_vel)
       IF (PRESENT(vel)) my_vel => vel
       mydt = 0.1_dp
       IF (PRESENT(dt)) mydt = dt
       CALL force_env_get(force_env,subsys=subsys,cell=cell,error=error)
       CALL cp_subsys_get(subsys, &
            atomic_kinds=atomic_kinds,&
            local_molecules_new=local_molecules,&
            local_particles=local_particles,&
            molecules_new=molecules,&
            molecule_kinds_new=molecule_kinds,&
            particles=particles,&
            gci=gci,&
            error=error)
       nparticle_kind = atomic_kinds%n_els
       has_vel=.FALSE.
       IF (.NOT.ASSOCIATED(my_vel)) THEN
          has_vel=.TRUE.
          ALLOCATE(my_vel(3,particles%n_els),stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          my_vel = 0.0_dp
          DO iparticle_kind=1,nparticle_kind
             nparticle_local = local_particles%n_el(iparticle_kind)
             DO iparticle_local=1,nparticle_local
                iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                my_vel (:,iparticle) = particles%els(iparticle)%v(:)
             END DO
          END DO
       END IF

       CALL rattle_control( gci=gci, local_molecules=local_molecules,&
            molecule_set=molecules%els, molecule_kind_set=molecule_kinds%els, &
            particle_set=particles%els, vel=my_vel, dt=mydt,&
            rattle_tol=shake_tol, log_unit=my_log_unit, lagrange_mult=my_lagrange_mult,&
            dump_lm=my_dump_lm, cell=cell, group=force_env%para_env%group,&
            local_particles=local_particles, error=error )

       ! Possibly reset the lagrange multipliers
       IF (PRESENT(reset)) THEN
          IF (reset) THEN
             ! Reset Intramolecular constraints
             DO i = 1,SIZE(molecules%els)
                IF (ASSOCIATED(molecules%els(i)%lci%lcolv)) THEN
                   DO j = 1, SIZE(molecules%els(i)%lci%lcolv)
                      ! Reset langrange multiplier
                      molecules%els(i)%lci%lcolv(j)%lambda = 0.0_dp
                   END DO
                END IF
                IF (ASSOCIATED(molecules%els(i)%lci%lg3x3)) THEN
                   DO j = 1, SIZE(molecules%els(i)%lci%lg3x3)
                      ! Reset langrange multiplier
                      molecules%els(i)%lci%lg3x3(j)%lambda = 0.0_dp
                   END DO
                END IF
                IF (ASSOCIATED(molecules%els(i)%lci%lg4x6)) THEN
                   DO j = 1, SIZE(molecules%els(i)%lci%lg4x6)
                      ! Reset langrange multiplier
                      molecules%els(i)%lci%lg4x6(j)%lambda = 0.0_dp
                   END DO
                END IF
             END DO
             ! Reset Intermolecular constraints
             IF (ASSOCIATED(gci)) THEN
                IF (ASSOCIATED(gci%lcolv)) THEN
                   DO j = 1, SIZE(gci%lcolv)
                      ! Reset langrange multiplier
                      gci%lcolv(j)%lambda = 0.0_dp
                   END DO
                END IF
                IF (ASSOCIATED(gci%lg3x3)) THEN
                   DO j = 1, SIZE(gci%lg3x3)
                      ! Reset langrange multiplier
                      gci%lg3x3(j)%lambda = 0.0_dp
                   END DO
                END IF
                IF (ASSOCIATED(gci%lg4x6)) THEN
                   DO j = 1, SIZE(gci%lg4x6)
                      ! Reset langrange multiplier
                      gci%lg4x6(j)%lambda = 0.0_dp
                   END DO
                END IF
             END IF
          END IF
       END IF

       IF (has_vel) THEN
          CALL update_particle_set ( particles%els, force_env%para_env%group, vel=my_vel)
       END IF
       DEALLOCATE(my_vel,stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE force_env_rattle
!***************************************************************************

!!****f* force_env_types/rescale_forces *
!!
!!   NAME
!!     rescale_forces
!!
!!   FUNCTION
!!     Rescale forces if requested
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - force_env: the force env to shake
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     tlaino
!!
!!*** **********************************************************************
  SUBROUTINE rescale_forces (force_env, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rescale_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iparticle
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp)                            :: force(3), max_value, mod_force
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: rescale_force_section

    failure=.FALSE.
    CALL timeset(routineN,handle)
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
    rescale_force_section => section_vals_get_subs_vals(force_env%force_env_section,"RESCALE_FORCES",error=error)
    CALL section_vals_get(rescale_force_section, explicit=explicit, error=error)
    IF (.NOT.failure.AND.explicit) THEN
       CALL section_vals_val_get(rescale_force_section,"MAX_FORCE",r_val=max_value,error=error)
       CALL force_env_get(force_env,subsys=subsys,error=error)
       CALL cp_subsys_get(subsys,particles=particles,error=error) 
       DO iparticle = 1, SIZE(particles%els)
          force = particles%els(iparticle)%f(:)
          mod_force = SQRT(DOT_PRODUCT(force,force))
          IF ((mod_force > max_value).AND.(mod_force /= 0.0_dp)) THEN
             force = force / mod_force * max_value 
             particles%els(iparticle)%f(:) = force
          END IF
       END DO
    END IF
  CALL timestop(handle)
END SUBROUTINE rescale_forces

END MODULE force_env_methods
