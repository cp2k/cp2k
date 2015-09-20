!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Interface for the force calculations
!> \par History
!>      cjm, FEB-20-2001: pass variable box_ref
!>      cjm, SEPT-12-2002: major reorganization
!>      fawzi, APR-12-2003: introduced force_env (based on the work by CJM&JGH)
!>      fawzi, NOV-3-2004: reorganized interface for f77 interface
!> \author fawzi
! *****************************************************************************
MODULE force_env_methods
  USE atprop_types,                    ONLY: atprop_init,&
                                             atprop_type
  USE cell_types,                      ONLY: cell_clone,&
                                             cell_create,&
                                             cell_release,&
                                             cell_type,&
                                             init_cell,&
                                             real_to_scaled,&
                                             scaled_to_real
  USE constraint_fxd,                  ONLY: fix_atom_control
  USE constraint_vsite,                ONLY: vsite_force_control
  USE cp_iter_types,                   ONLY: cp_iteration_info_copy_iter
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr,&
                                             low_print_level
  USE cp_para_env,                     ONLY: cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_result_methods,               ONLY: cp_results_erase,&
                                             cp_results_mp_bcast,&
                                             get_results,&
                                             test_for_result
  USE cp_result_types,                 ONLY: cp_result_copy,&
                                             cp_result_create,&
                                             cp_result_p_type,&
                                             cp_result_release,&
                                             cp_result_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_p_type,&
                                             cp_subsys_set,&
                                             cp_subsys_type
  USE eip_environment_types,           ONLY: eip_env_retain,&
                                             eip_environment_type
  USE eip_silicon,                     ONLY: eip_bazant,&
                                             eip_lenosky
  USE external_potential_methods,      ONLY: add_external_potential
  USE fist_environment_types,          ONLY: fist_env_retain,&
                                             fist_environment_type
  USE fist_force,                      ONLY: fist_calc_energy_force
  USE force_env_types,                 ONLY: &
       force_env_get, force_env_get_natom, force_env_p_type, force_env_set, &
       force_env_type, use_eip_force, use_fist_force, use_mixed_force, &
       use_prog_name, use_qmmm, use_qmmmx, use_qs_force
  USE force_env_utils,                 ONLY: rescale_forces,&
                                             write_forces,&
                                             write_stress_tensor
  USE force_fields_util,               ONLY: get_generic_info
  USE fp_methods,                      ONLY: fp_eval
  USE fparser,                         ONLY: EvalErrType,&
                                             evalf,&
                                             evalfd,&
                                             finalizef,&
                                             initf,&
                                             parsef
  USE global_types,                    ONLY: global_environment_type,&
                                             globenv_retain
  USE input_constants,                 ONLY: debug_run,&
                                             mix_coupled,&
                                             mix_generic,&
                                             mix_linear_combination,&
                                             mix_minimum,&
                                             mix_restrained,&
                                             use_bazant_eip,&
                                             use_lenosky_eip
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_retain,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE machine,                         ONLY: m_memory
  USE message_passing,                 ONLY: mp_sum,&
                                             mp_sync
  USE metadynamics_types,              ONLY: meta_env_retain,&
                                             meta_env_type
  USE mixed_energy_types,              ONLY: mixed_energy_type,&
                                             mixed_force_type
  USE mixed_environment_types,         ONLY: get_mixed_env,&
                                             mixed_env_retain,&
                                             mixed_environment_type
  USE mixed_environment_utils,         ONLY: get_subsys_map_index,&
                                             mixed_map_forces
  USE particle_list_types,             ONLY: particle_list_p_type,&
                                             particle_list_type
  USE physcon,                         ONLY: debye
  USE qmmm_force,                      ONLY: qmmm_calc_energy_force
  USE qmmm_types,                      ONLY: qmmm_env_retain,&
                                             qmmm_env_type
  USE qmmmx_force,                     ONLY: qmmmx_calc_energy_force
  USE qmmmx_types,                     ONLY: qmmmx_env_retain,&
                                             qmmmx_env_type
  USE qs_environment_types,            ONLY: qs_env_retain,&
                                             qs_environment_type
  USE qs_force,                        ONLY: qs_calc_energy_force
  USE restraint,                       ONLY: restraint_control
  USE string_utilities,                ONLY: compress
  USE termination,                     ONLY: print_warning
  USE virial_types,                    ONLY: cp_virial,&
                                             sym_virial,&
                                             virial_create,&
                                             virial_p_type,&
                                             virial_release,&
                                             virial_type,&
                                             zero_virial
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'force_env_methods'

  PUBLIC :: force_env_create,&
            force_env_calc_energy_force,&
            force_env_calc_num_pressure

  INTEGER, SAVE, PRIVATE :: last_force_env_id=0

CONTAINS

! *****************************************************************************
!> \brief Interface routine for force and energy calculations
!> \param force_env the force_env of which you want the energy and forces
!> \param calc_force if false the forces *might* be left unchanged
!>        or be unvalid, no guarantee on them is done.Defaults to true
!> \param consistent_energies Performs an additional qs_ks_update_qs_env, so
!>          that the energies are appropriate to the forces, they are in the
!>          non-selfconsistent case not consistent to each other! [08.2005, TdK]
!> \param skip_external_control ...
!> \param eval_energy_forces ...
!> \param require_consistent_energy_force ...
!> \param linres ...
!> \param calc_stress_tensor ...
!> \author CJM & fawzi
! *****************************************************************************
  RECURSIVE SUBROUTINE force_env_calc_energy_force ( force_env, calc_force, &
       consistent_energies, skip_external_control, eval_energy_forces, &
       require_consistent_energy_force, linres, calc_stress_tensor)

    TYPE(force_env_type), POINTER            :: force_env
    LOGICAL, INTENT(IN), OPTIONAL :: calc_force, consistent_energies, &
      skip_external_control, eval_energy_forces, &
      require_consistent_energy_force, linres, calc_stress_tensor

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_calc_energy_force', &
      routineP = moduleN//':'//routineN
    REAL(kind=dp), PARAMETER                 :: ateps = 1.0E-6_dp

    INTEGER                                  :: i, j, nat, ndigits, &
                                                output_unit, print_forces
    LOGICAL :: calculate_forces, calculate_stress_tensor, energy_consistency, &
      eval_ef, failure, linres_run, my_skip
    REAL(KIND=dp)                            :: checksum, e_pot, sum_energy, &
                                                sum_pv_virial, &
                                                sum_stress_tensor
    REAL(KIND=dp), DIMENSION(3)              :: grand_total_force, total_force
    REAL(KIND=dp), DIMENSION(3, 3)           :: atomic_stress_tensor, &
                                                diff_stress_tensor
    TYPE(atprop_type), POINTER               :: atprop_env
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(virial_type), POINTER               :: virial

    NULLIFY (logger, virial, subsys, atprop_env, cell)
    logger => cp_get_default_logger()
    failure = .FALSE.
    eval_ef = .TRUE.
    my_skip = .FALSE.
    calculate_forces = .TRUE.
    energy_consistency = .FALSE.
    linres_run = .FALSE.
    IF (PRESENT(eval_energy_forces)) eval_ef = eval_energy_forces
    IF (PRESENT(skip_external_control)) my_skip = skip_external_control
    IF (PRESENT(calc_force)) calculate_forces = calc_force
    IF (PRESENT(calc_stress_tensor)) THEN
       calculate_stress_tensor = calc_stress_tensor
    ELSE
       calculate_stress_tensor = calculate_forces
    END IF
    IF (PRESENT(consistent_energies)) energy_consistency = consistent_energies
    IF (PRESENT(linres)) linres_run = linres

    CPASSERT(ASSOCIATED(force_env))
    CPASSERT(force_env%ref_count>0)
    CALL force_env_get(force_env,subsys=subsys)
    CALL force_env_set(force_env,additional_potential=0.0_dp)
    CALL cp_subsys_get(subsys, virial=virial, atprop=atprop_env, cell=cell)
    IF (virial%pv_availability) CALL zero_virial(virial,reset=.FALSE.)

    nat=force_env_get_natom(force_env)
    CALL atprop_init(atprop_env,nat)
    IF (eval_ef) THEN
       SELECT CASE ( force_env%in_use )
       CASE ( use_fist_force )
          CALL fist_calc_energy_force(force_env%fist_env)
       CASE ( use_qs_force )
          CALL qs_calc_energy_force(force_env%qs_env,calculate_forces,energy_consistency,linres_run)
       CASE (use_eip_force)
          IF (force_env%eip_env%eip_model == use_lenosky_eip) THEN
             CALL eip_lenosky(force_env%eip_env)
          ELSE IF (force_env%eip_env%eip_model == use_bazant_eip) THEN
             CALL eip_bazant(force_env%eip_env)
          END IF
       CASE ( use_qmmm )
          CALL qmmm_calc_energy_force(force_env%qmmm_env,&
               calculate_forces,energy_consistency,linres=linres_run)
       CASE ( use_qmmmx )
          CALL qmmmx_calc_energy_force(force_env%qmmmx_env,&
               calculate_forces,energy_consistency,linres=linres_run,&
               require_consistent_energy_force=require_consistent_energy_force)
       CASE ( use_mixed_force )
          CALL mixed_energy_forces(force_env,calculate_forces)
       CASE default
          CPASSERT(.FALSE.)
       END SELECT
    END IF
    ! In case it is requested, we evaluate the stress tensor numerically
    IF (virial%pv_availability) THEN
       IF (virial%pv_numer.AND.calculate_stress_tensor) THEN
          ! Compute the numerical stress tensor
          CALL force_env_calc_num_pressure(force_env)
       ELSE
          IF (calculate_forces) THEN
             ! Symmetrize analytical stress tensor
             CALL sym_virial(virial)
          ELSE
             IF (calculate_stress_tensor) THEN
                CALL print_warning(routineN,moduleN,__LINE__,&
                                   "The calculation of the stress tensor "//&
                                   "requires the calculation of the forces",&
                                   force_env%para_env)
             END IF
          END IF
       END IF
    END IF

    !sample peak memory
    CALL m_memory()

    ! Some additional tasks..
    IF (.NOT.my_skip) THEN
       ! Flexible Partitioning
       IF (ASSOCIATED(force_env%fp_env)) THEN
          IF (force_env%fp_env%use_fp) THEN
             CALL fp_eval(force_env%fp_env,subsys,cell)
          ENDIF
       ENDIF
       ! Constraints ONLY of Fixed Atom type
       CALL fix_atom_control(force_env)
       ! All Restraints
       CALL restraint_control(force_env)
       ! Virtual Sites
       CALL vsite_force_control(force_env)
       ! External Potential
       CALL add_external_potential(force_env)
       ! Rescale forces if requested
       CALL rescale_forces(force_env)
    END IF

    ! Print always Energy in the same format for all methods
    output_unit = cp_print_key_unit_nr(logger,force_env%force_env_section,"PRINT%PROGRAM_RUN_INFO",&
      extension=".Log")
    IF (output_unit > 0) THEN
       CALL force_env_get(force_env, potential_energy=e_pot)
       WRITE(output_unit,'(/,T2,"ENERGY| Total FORCE_EVAL ( ",A," ) energy (a.u.): ",T55,F26.15,/)')&
            ADJUSTR(TRIM(use_prog_name(force_env%in_use))),e_pot
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,force_env%force_env_section,&
                                      "PRINT%PROGRAM_RUN_INFO")

    ! Print forces, if requested
    print_forces = cp_print_key_unit_nr(logger,force_env%force_env_section,"PRINT%FORCES",&
                                        extension=".xyz")
    IF ((print_forces > 0).AND.calculate_forces) THEN
       CALL force_env_get(force_env,subsys=subsys)
       CALL cp_subsys_get(subsys,&
                          core_particles=core_particles,&
                          particles=particles,&
                          shell_particles=shell_particles)
       ! Variable precision output of the forces
       CALL section_vals_val_get(force_env%force_env_section,"PRINT%FORCES%NDIGITS",&
                                 i_val=ndigits)
       CALL write_forces(particles,print_forces,"ATOMIC",ndigits,total_force)
       grand_total_force(:) = total_force(:)
       IF (ASSOCIATED(core_particles)) THEN
          CALL write_forces(core_particles,print_forces,"CORE",ndigits,total_force)
          grand_total_force(:) = grand_total_force(:) + total_force(:)
       END IF
       IF (ASSOCIATED(shell_particles)) THEN
          CALL write_forces(shell_particles,print_forces,"SHELL",ndigits,total_force,&
                            grand_total_force)
       END IF
    END IF
    CALL cp_print_key_finished_output(print_forces,logger,force_env%force_env_section,"PRINT%FORCES")


    ! Write stress tensor
    IF (virial%pv_availability) THEN
       ! If the virial is defined but we are not computing forces let's zero the
       ! virial for consistency
       IF (calculate_forces.AND.calculate_stress_tensor) THEN
          output_unit = cp_print_key_unit_nr(logger,force_env%force_env_section,"PRINT%STRESS_TENSOR",&
                                             extension=".stress_tensor")
          IF (output_unit > 0) THEN
             CALL section_vals_val_get(force_env%force_env_section,"PRINT%STRESS_TENSOR%NDIGITS",&
                                       i_val=ndigits)
             CALL write_stress_tensor(virial%pv_virial,output_unit,cell,ndigits,virial%pv_numer)
          END IF
          CALL cp_print_key_finished_output(output_unit,logger,force_env%force_env_section,&
                                            "PRINT%STRESS_TENSOR")
       ELSE
          CALL zero_virial(virial,reset=.FALSE.)
       END IF
    ELSE
       output_unit = cp_print_key_unit_nr(logger,force_env%force_env_section,"PRINT%STRESS_TENSOR",&
                                          extension=".stress_tensor")
       IF (output_unit > 0) THEN
          CALL print_warning(routineN,moduleN,__LINE__,"To print the stress tensor switch on the "//&
                             "virial evaluation with the keyword: STRESS_TENSOR",force_env%para_env)
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,force_env%force_env_section,&
                                         "PRINT%STRESS_TENSOR")
    END IF

    ! Atomic energy      
    output_unit = cp_print_key_unit_nr(logger,force_env%force_env_section,"PRINT%PROGRAM_RUN_INFO",&
                                       extension=".Log")
    IF (atprop_env%energy) THEN
       CALL mp_sum(atprop_env%atener, force_env%para_env%group)
       CALL force_env_get(force_env, potential_energy=e_pot)
       IF (output_unit > 0) THEN
          IF (logger%iter_info%print_level > low_print_level) THEN
             WRITE (UNIT=output_unit,FMT="(/,T6,A,T15,A)") "Atom","Potential energy"
             WRITE (UNIT=output_unit,FMT="(T2,I8,1X,F20.10)")&
              (i,atprop_env%atener(i),i=1,SIZE(atprop_env%atener))
          END IF
          sum_energy = accurate_sum(atprop_env%atener(:))
          checksum = ABS(e_pot - sum_energy)
          WRITE (UNIT=output_unit,FMT="(/,(T2,A,T56,F25.13))")&
           "Potential energy (Atomic):",sum_energy,&
           "Potential energy (Total) :",e_pot,&
           "Difference               :",checksum
          CPASSERT((checksum < ateps*ABS(e_pot)))
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,force_env%force_env_section,&
                                         "PRINT%PROGRAM_RUN_INFO")
    END IF

    ! Atomic stress    
    output_unit = cp_print_key_unit_nr(logger,force_env%force_env_section,"PRINT%PROGRAM_RUN_INFO",&
                                       extension=".Log")

    IF (atprop_env%stress) THEN
       CALL mp_sum(atprop_env%atstress,force_env%para_env%group)
       ! symmetrize (same as pv_virial)
       DO i=1,SIZE(atprop_env%atstress,3)
          atprop_env%atstress(:,:,i) = 0.5_dp * (atprop_env%atstress(:,:,i)  &
                                              + TRANSPOSE(atprop_env%atstress(:,:,i)))
       END DO
       IF (output_unit > 0) THEN
          IF (logger%iter_info%print_level > low_print_level) THEN
             DO i=1,SIZE(atprop_env%atstress,3)
                WRITE (UNIT=output_unit,FMT="(/,T2,I0,T16,A1,2(19X,A1))") i,"X","Y","Z"
                WRITE (UNIT=output_unit,FMT="(A3,3F20.13)") "X",(atprop_env%atstress(1,j,i),j=1,3)
                WRITE (UNIT=output_unit,FMT="(A3,3F20.13)") "Y",(atprop_env%atstress(2,j,i),j=1,3)
                WRITE (UNIT=output_unit,FMT="(A3,3F20.13)") "Z",(atprop_env%atstress(3,j,i),j=1,3)
                WRITE (UNIT=output_unit,FMT="(T2,A,F20.13)") "1/3 Trace(Atomic stress tensor):",&
                 (atprop_env%atstress(1,1,i) + atprop_env%atstress(2,2,i) + atprop_env%atstress(3,3,i))/3.0_dp
             END DO
          END IF
          atomic_stress_tensor(:,:) = 0.0_dp
          DO i=1,3
             atomic_stress_tensor(i,i) = accurate_sum(atprop_env%atstress(i,i,:))
             DO j=i+1,3
                atomic_stress_tensor(i,j) = accurate_sum(atprop_env%atstress(i,j,:))
                atomic_stress_tensor(j,i) = atomic_stress_tensor(i,j)
             END DO
          END DO
          WRITE (UNIT=output_unit,FMT="(/,T2,A,T15,A1,2(19X,A1))") "Atomic","X","Y","Z"
          WRITE (UNIT=output_unit,FMT="(A3,3F20.13)") "X",(atomic_stress_tensor(1,i),i=1,3)
          WRITE (UNIT=output_unit,FMT="(A3,3F20.13)") "Y",(atomic_stress_tensor(2,i),i=1,3)
          WRITE (UNIT=output_unit,FMT="(A3,3F20.13)") "Z",(atomic_stress_tensor(3,i),i=1,3)
          WRITE (UNIT=output_unit,FMT="(T2,A,10X,F20.13)") "1/3 Trace(Atomic stress tensor):",&
           (atomic_stress_tensor(1,1) + atomic_stress_tensor(2,2) + atomic_stress_tensor(3,3))/3.0_dp
          sum_stress_tensor = accurate_sum(atomic_stress_tensor(:,:))
          IF (virial%pv_availability.AND.calculate_forces) THEN
             WRITE (UNIT=output_unit,FMT="(/,T2,A,T16,A1,2(19X,A1))") "Total","X","Y","Z"
             WRITE (UNIT=output_unit,FMT="(A3,3F20.13)") "X",(virial%pv_virial(1,i),i=1,3)
             WRITE (UNIT=output_unit,FMT="(A3,3F20.13)") "Y",(virial%pv_virial(2,i),i=1,3)
             WRITE (UNIT=output_unit,FMT="(A3,3F20.13)") "Z",(virial%pv_virial(3,i),i=1,3)
             WRITE (UNIT=output_unit,FMT="(T2,A,10X,F20.13)") "1/3 Trace(Total stress tensor): ",&
              (virial%pv_virial(1,1) + virial%pv_virial(2,2) + virial%pv_virial(3,3))/3.0_dp
             sum_pv_virial = SUM(virial%pv_virial(:,:))
             diff_stress_tensor(:,:) = ABS(virial%pv_virial(:,:) - atomic_stress_tensor(:,:))
             WRITE (UNIT=output_unit,FMT="(/,T2,A,T16,A1,2(19X,A1))") "Diff","X","Y","Z"
             WRITE (UNIT=output_unit,FMT="(A3,3F20.13)") "X",(diff_stress_tensor(1,i),i=1,3)
             WRITE (UNIT=output_unit,FMT="(A3,3F20.13)") "Y",(diff_stress_tensor(2,i),i=1,3)
             WRITE (UNIT=output_unit,FMT="(A3,3F20.13)") "Z",(diff_stress_tensor(3,i),i=1,3)
             WRITE (UNIT=output_unit,FMT="(T2,A,10X,F20.13)") "1/3 Trace(Diff)               : ",&
              (diff_stress_tensor(1,1) + diff_stress_tensor(2,2) + diff_stress_tensor(3,3))/3.0_dp
             checksum = accurate_sum(diff_stress_tensor(:,:))
             WRITE (UNIT=output_unit,FMT="(/,(T2,A,11X,F25.13))")&
              "Checksum stress (Atomic) :",sum_stress_tensor,&
              "Checksum stress (Total)  :",sum_pv_virial,&
              "Difference               :",checksum
             CPASSERT(checksum < ateps)
          END IF
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,force_env%force_env_section,&
                                         "PRINT%PROGRAM_RUN_INFO")
    END IF

  END SUBROUTINE force_env_calc_energy_force

! *****************************************************************************
!> \brief Evaluates the stress tensor and pressure numerically
!> \param force_env ...
!> \param dx ...
!> \par History
!>      10.2005 created [JCS]
!>      05.2009 Teodoro Laino [tlaino] - rewriting for general force_env
!>
!> \author JCS
! *****************************************************************************
  SUBROUTINE force_env_calc_num_pressure(force_env,dx)

    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: dx

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_calc_num_pressure', &
      routineP = moduleN//':'//routineN
    REAL(kind=dp), PARAMETER                 :: default_dx = 0.001_dp

    INTEGER                                  :: i, ip, iq, istat, j, k, &
                                                natom, ncore, ndigits, &
                                                nshell, output_unit
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dx_w
    REAL(KIND=dp), DIMENSION(2)              :: numer_energy
    REAL(KIND=dp), DIMENSION(3)              :: s
    REAL(KIND=dp), DIMENSION(3, 3)           :: numer_stress
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ref_pos_atom, ref_pos_core, &
                                                ref_pos_shell
    TYPE(cell_type), POINTER                 :: cell, cell_local
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(virial_type), POINTER               :: virial

    NULLIFY (cell_local)
    NULLIFY (core_particles)
    NULLIFY (particles)
    NULLIFY (shell_particles)
    NULLIFY (ref_pos_atom)
    NULLIFY (ref_pos_core)
    NULLIFY (ref_pos_shell)
    natom = 0
    ncore = 0
    nshell = 0
    numer_stress = 0.0_dp

    failure = .FALSE.
    logger => cp_get_default_logger()


    dx_w = default_dx
    IF (PRESENT(dx)) dx_w = dx
    CALL force_env_get(force_env,subsys=subsys,globenv=globenv)
    CALL cp_subsys_get(subsys,&
                       core_particles=core_particles,&
                       particles=particles,&
                       shell_particles=shell_particles,&
                       virial=virial)
    output_unit = cp_print_key_unit_nr(logger,force_env%force_env_section,"PRINT%STRESS_TENSOR",&
                                       extension=".stress_tensor")
    IF (output_unit > 0) THEN
       WRITE (output_unit,'(/A,A/)') ' **************************** ', &
         'NUMERICAL STRESS ********************************'
    END IF

    ! Save all original particle positions
    natom = particles%n_els
    ALLOCATE (ref_pos_atom(natom,3),STAT=istat)
    CPASSERT(istat == 0)
    DO i=1,natom
       ref_pos_atom(i,:) = particles%els(i)%r
    END DO
    IF (ASSOCIATED(core_particles)) THEN
       ncore = core_particles%n_els
       ALLOCATE (ref_pos_core(ncore,3),STAT=istat)
       CPASSERT(istat == 0)
       DO i=1,ncore
          ref_pos_core(i,:) = core_particles%els(i)%r
       END DO
    END IF
    IF (ASSOCIATED(shell_particles)) THEN
       nshell = shell_particles%n_els
       ALLOCATE (ref_pos_shell(nshell,3),STAT=istat)
       CPASSERT(istat == 0)
       DO i=1,nshell
          ref_pos_shell(i,:) = shell_particles%els(i)%r
       END DO
    END IF
    CALL force_env_get(force_env,cell=cell)
    CALL cell_create(cell_local)
    CALL cell_clone(cell,cell_local)
    ! First change box
    DO ip=1,3
       DO iq=1,3
          IF (virial%pv_diagonal.AND.(ip /= iq)) CYCLE
          DO k=1,2
             cell%hmat(ip,iq) = cell_local%hmat(ip,iq) - (-1.0_dp)**k*dx_w
             CALL init_cell(cell)
             ! Scale positions
             DO i=1,natom
                CALL real_to_scaled(s,ref_pos_atom(i,1:3),cell_local)
                CALL scaled_to_real(particles%els(i)%r,s,cell)
             END DO
             DO i=1,ncore
                CALL real_to_scaled(s,ref_pos_core(i,1:3),cell_local)
                CALL scaled_to_real(core_particles%els(i)%r,s,cell)
             END DO
             DO i=1,nshell
                CALL real_to_scaled(s,ref_pos_shell(i,1:3),cell_local)
                CALL scaled_to_real(shell_particles%els(i)%r,s,cell)
             END DO

             ! Compute energies
             CALL force_env_calc_energy_force(force_env,&
                                              calc_force=.FALSE.,&
                                              consistent_energies=.TRUE.,&
                                              calc_stress_tensor=.FALSE.)
             CALL force_env_get(force_env,potential_energy=numer_energy(k))
             ! Reset cell
             cell%hmat(ip,iq) = cell_local%hmat(ip,iq)
          END DO
          CALL init_cell(cell)
          numer_stress(ip,iq) = 0.5_dp*(numer_energy(1) - numer_energy(2))/dx_w
          IF (output_unit > 0) THEN
             WRITE (UNIT=output_unit,FMT="(T7,A,F7.4,A,T27,A,F7.4,A,T49,A)")&
              "E("//ACHAR(119+ip)//ACHAR(119+iq)//" +",dx_w,")",&
              "E("//ACHAR(119+ip)//ACHAR(119+iq)//" -",dx_w,")",&
              "f(numerical)"
             WRITE (UNIT=output_unit,FMT="(3(1X,F19.8))")&
              numer_energy(1:2),numer_stress(ip,iq)
          END IF
       END DO
    END DO

    ! Reset positions and rebuild original environment
    DO i=1,natom
       particles%els(i)%r = ref_pos_atom(i,:)
    END DO
    DO i=1,ncore
       core_particles%els(i)%r = ref_pos_core(i,:)
    END DO
    DO i=1,nshell
       shell_particles%els(i)%r = ref_pos_shell(i,:)
    END DO
    CALL force_env_calc_energy_force(force_env,&
                                     calc_force=.FALSE.,&
                                     consistent_energies=.TRUE.,&
                                     calc_stress_tensor=.FALSE.)

    ! Computing pv_test
    virial%pv_virial = 0.0_dp
    DO i=1,3
       DO j=1,3
          DO k=1,3
             virial%pv_virial(i,j) = virial%pv_virial(i,j) -&
                                     0.5_dp*(numer_stress(i,k)*cell_local%hmat(j,k) +&
                                             numer_stress(j,k)*cell_local%hmat(i,k))
          END DO
       END DO
    END DO

    IF (output_unit > 0) THEN
       IF (globenv%run_type_id == debug_run) THEN
          CALL section_vals_val_get(force_env%force_env_section,"PRINT%STRESS_TENSOR%NDIGITS",&
                                    i_val=ndigits)
          CALL write_stress_tensor(virial%pv_virial,output_unit,cell,ndigits,virial%pv_numer)
       END IF
       WRITE (output_unit,'(/,A,/)') ' **************************** '//&
         'NUMERICAL STRESS END *****************************'
    END IF

    CALL cp_print_key_finished_output(output_unit,logger,force_env%force_env_section,&
                                      "PRINT%STRESS_TENSOR")

    ! Release storage
    IF (ASSOCIATED(ref_pos_atom)) THEN
       DEALLOCATE (ref_pos_atom,STAT=istat)
       CPASSERT(istat == 0)
    END IF
    IF (ASSOCIATED(ref_pos_core)) THEN
       DEALLOCATE (ref_pos_core,STAT=istat)
       CPASSERT(istat == 0)
    END IF
    IF (ASSOCIATED(ref_pos_shell)) THEN
       DEALLOCATE (ref_pos_shell,STAT=istat)
       CPASSERT(istat == 0)
    END IF
    IF (ASSOCIATED(cell_local)) CALL cell_release(cell_local)


  END SUBROUTINE force_env_calc_num_pressure

! *****************************************************************************
!> \brief creates and initializes a force environment
!> \param force_env the force env to create
!> \param root_section ...
!> \param para_env ...
!> \param globenv ...
!> \param fist_env , qs_env: exactly one of these should be
!>        associated, the one that is active
!> \param qs_env ...
!> \param meta_env ...
!> \param sub_force_env ...
!> \param qmmm_env ...
!> \param qmmmx_env ...
!> \param eip_env ...
!> \param force_env_section ...
!> \param mixed_env ...
!> \par History
!>      04.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE force_env_create(force_env,root_section,para_env,globenv,fist_env,&
       qs_env,meta_env,sub_force_env,qmmm_env,qmmmx_env,eip_env,force_env_section,&
       mixed_env)

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(fist_environment_type), OPTIONAL, &
      POINTER                                :: fist_env
    TYPE(qs_environment_type), OPTIONAL, &
      POINTER                                :: qs_env
    TYPE(meta_env_type), OPTIONAL, POINTER   :: meta_env
    TYPE(force_env_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: sub_force_env
    TYPE(qmmm_env_type), OPTIONAL, POINTER   :: qmmm_env
    TYPE(qmmmx_env_type), OPTIONAL, POINTER  :: qmmmx_env
    TYPE(eip_environment_type), OPTIONAL, &
      POINTER                                :: eip_env
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(mixed_environment_type), OPTIONAL, &
      POINTER                                :: mixed_env

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE ( force_env, stat=stat )
    CPASSERT(stat==0)
    NULLIFY (force_env%fist_env, force_env%qs_env, &
         force_env%para_env, force_env%globenv, &
         force_env%meta_env, force_env%sub_force_env, &
         force_env%qmmm_env, force_env%qmmmx_env, force_env%fp_env, &
         force_env%force_env_section, force_env%eip_env,force_env%mixed_env,&
         force_env%root_section)
    last_force_env_id=last_force_env_id+1
    force_env%id_nr=last_force_env_id
    force_env%ref_count=1
    force_env%in_use=0
    force_env%additional_potential=0.0_dp

    force_env%globenv => globenv
    CALL globenv_retain(force_env%globenv)

    force_env%root_section => root_section
    CALL section_vals_retain(root_section)

    force_env%para_env=>para_env
    CALL cp_para_env_retain(force_env%para_env)

    CALL section_vals_retain(force_env_section)
    force_env%force_env_section => force_env_section

    IF (PRESENT(fist_env)) THEN
       CPASSERT(ASSOCIATED(fist_env))
       CPASSERT(force_env%in_use==0)
       force_env%in_use=use_fist_force
       force_env%fist_env => fist_env
       CALL fist_env_retain(fist_env)
    END IF
    IF (PRESENT(eip_env)) THEN
       CPASSERT(ASSOCIATED(eip_env))
       CPASSERT(force_env%in_use==0)
       force_env%in_use = use_eip_force
       force_env%eip_env => eip_env
       CALL eip_env_retain(eip_env)
    END IF
    IF (PRESENT(qs_env)) THEN
       CPASSERT(ASSOCIATED(qs_env))
       CPASSERT(force_env%in_use==0)
       force_env%in_use=use_qs_force
       force_env%qs_env => qs_env
       CALL qs_env_retain(qs_env)
    END IF
    IF (PRESENT(qmmm_env)) THEN
       CPASSERT(ASSOCIATED(qmmm_env))
       CPASSERT(force_env%in_use==0)
       force_env%in_use=use_qmmm
       force_env%qmmm_env => qmmm_env
       CALL qmmm_env_retain(qmmm_env)
    END IF
    IF (PRESENT(qmmmx_env)) THEN
       CPASSERT(ASSOCIATED(qmmmx_env))
       CPASSERT(force_env%in_use==0)
       force_env%in_use=use_qmmmx
       force_env%qmmmx_env => qmmmx_env
       CALL qmmmx_env_retain(qmmmx_env)
    END IF
    IF (PRESENT(mixed_env)) THEN
       CPASSERT(ASSOCIATED(mixed_env))
       CPASSERT(force_env%in_use==0)
       force_env%in_use=use_mixed_force
       force_env%mixed_env => mixed_env
       CALL mixed_env_retain ( mixed_env)
    END IF
    CPASSERT(force_env%in_use/=0)

    IF (PRESENT(sub_force_env)) THEN
       force_env%sub_force_env => sub_force_env
    END IF

    IF (PRESENT(meta_env)) THEN
       force_env%meta_env => meta_env
       CALL meta_env_retain(meta_env)
    ELSE
       NULLIFY(force_env%meta_env)
    END IF


  END SUBROUTINE force_env_create


! *****************************************************************************
!> \brief ****f* force_env_methods/mixed_energy_forces  [1.0]
!>
!>     Computes energy and forces for a mixed force_env type
!> \param force_env ...
!> \param calculate_forces ...
!> \par History
!>       11.06  created [fschiff]
!>       04.07  generalization to an illimited number of force_eval [tlaino]
!>       04.07  further generalization to force_eval with different geometrical
!>              structures [tlaino]
!>       04.08  reorganizing the genmix structure (collecting common code)
!> \author Florian Schiffmann
! *****************************************************************************
  SUBROUTINE mixed_energy_forces(force_env, calculate_forces)

    TYPE(force_env_type), POINTER            :: force_env
    LOGICAL, INTENT(IN)                      :: calculate_forces

    CHARACTER(len=*), PARAMETER :: routineN = 'mixed_energy_forces', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: coupling_function
    CHARACTER(LEN=default_string_length)     :: def_error, description, &
                                                this_error
    INTEGER :: iforce_eval, iparticle, jparticle, mixing_type, my_group, &
      natom, nforce_eval, source, stat, unit_nr
    INTEGER, DIMENSION(:), POINTER           :: glob_natoms, map_index
    LOGICAL                                  :: dip_exists, failure
    REAL(KIND=dp) :: coupling_parameter, dedf, der_1, der_2, dx, energy, err, &
      lambda, lerr, restraint_strength, restraint_target, sd
    REAL(KIND=dp), DIMENSION(3)              :: dip_mix
    REAL(KIND=dp), DIMENSION(:), POINTER     :: energies
    TYPE(cell_type), POINTER                 :: cell_mix
    TYPE(cp_logger_type), POINTER            :: logger, my_logger
    TYPE(cp_result_p_type), DIMENSION(:), &
      POINTER                                :: results
    TYPE(cp_result_type), POINTER            :: loc_results, results_mix
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
                                                gen_section, mapping_section, &
                                                mixed_section, root_section
    TYPE(virial_p_type), DIMENSION(:), &
      POINTER                                :: virials
    TYPE(virial_type), POINTER               :: loc_virial, virial_mix

    failure=.FALSE.
    logger => cp_get_default_logger()
    CPASSERT(ASSOCIATED(force_env))
    ! Get infos about the mixed subsys
    CALL force_env_get(force_env=force_env,&
                       subsys=subsys_mix,&
                       force_env_section=force_env_section,&
                       root_section=root_section,&
                       cell=cell_mix)
    CALL cp_subsys_get(subsys=subsys_mix,&
                       particles=particles_mix,&
                       virial=virial_mix,&
                       results=results_mix)
    NULLIFY(map_index, glob_natoms, global_forces)

    nforce_eval = SIZE(force_env%sub_force_env)
    mixed_section => section_vals_get_subs_vals(force_env_section,"MIXED")
    mapping_section => section_vals_get_subs_vals(mixed_section,"MAPPING")
    ! Global Info
    ALLOCATE(subsystems(nforce_eval), stat=stat)
    CPASSERT(stat==0)
    ALLOCATE(particles(nforce_eval), stat=stat)
    CPASSERT(stat==0)
    ! Local Info to sync
    ALLOCATE(global_forces(nforce_eval), stat=stat)
    CPASSERT(stat==0)
    ALLOCATE(energies(nforce_eval), stat=stat)
    CPASSERT(stat==0)
    ALLOCATE(glob_natoms(nforce_eval), stat=stat)
    CPASSERT(stat==0)
    ALLOCATE(virials(nforce_eval), stat=stat)
    CPASSERT(stat==0)
    ALLOCATE(results(nforce_eval), stat=stat)
    CPASSERT(stat==0)
    energies    = 0.0_dp
    glob_natoms = 0
    DO iforce_eval = 1, nforce_eval
       NULLIFY(subsystems(iforce_eval)%subsys, particles(iforce_eval)%list)
       NULLIFY(results(iforce_eval)%results, virials(iforce_eval)%virial)
       CALL virial_create (virials(iforce_eval)%virial)
       CALL cp_result_create (results(iforce_eval)%results)
       IF (.NOT.ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) CYCLE
       ! From this point on the error is the sub_error
       my_group = force_env%mixed_env%group_distribution(force_env%para_env%mepos)
       my_logger => force_env%mixed_env%sub_logger(my_group+1)%p
       ! Copy iterations info (they are updated only in the main mixed_env)
       CALL cp_iteration_info_copy_iter(logger%iter_info, my_logger%iter_info)
       CALL cp_add_default_logger(my_logger)

       ! Get all available subsys
       CALL force_env_get(force_env=force_env%sub_force_env(iforce_eval)%force_env,&
                          subsys=subsystems(iforce_eval)%subsys)

       ! all force_env share the same cell
       CALL cp_subsys_set(subsystems(iforce_eval)%subsys, cell=cell_mix)

       ! Get available particles
       CALL cp_subsys_get(subsys=subsystems(iforce_eval)%subsys,&
                          particles=particles(iforce_eval)%list)

       ! Get Mapping index array
       natom = SIZE(particles(iforce_eval)%list%els)
       CALL get_subsys_map_index(mapping_section, natom, iforce_eval, nforce_eval, &
            map_index)

       ! Mapping particles from iforce_eval environment to the mixed env
       DO iparticle = 1, natom
          jparticle = map_index(iparticle)
          particles(iforce_eval)%list%els(iparticle)%r= particles_mix%els(jparticle)%r
       END DO

       ! Calculate energy and forces for each sub_force_env
       CALL force_env_calc_energy_force(force_env%sub_force_env(iforce_eval)%force_env,&
                                        calc_force=calculate_forces,&
                                        skip_external_control=.TRUE.)
       ! Only the rank 0 process collect info for each computation
       IF ( force_env%sub_force_env(iforce_eval)%force_env%para_env%mepos==&
            force_env%sub_force_env(iforce_eval)%force_env%para_env%source) THEN
          CALL force_env_get(force_env%sub_force_env(iforce_eval)%force_env,&
                             potential_energy=energy)
          CALL cp_subsys_get(subsystems(iforce_eval)%subsys,&
                             virial=loc_virial, results=loc_results)
          energies(iforce_eval)    = energy
          glob_natoms(iforce_eval) = natom
          CALL cp_virial(loc_virial, virials(iforce_eval)%virial)
          CALL cp_result_copy(loc_results, results(iforce_eval)%results)
       END IF
       ! Deallocate map_index array
       IF (ASSOCIATED(map_index)) THEN
          DEALLOCATE(map_index, stat=stat)
          CPASSERT(stat==0)
       END IF
       CALL cp_rm_default_logger()
    END DO

    ! Handling Parallel execution
    CALL mp_sync(force_env%para_env%group)
    ! Let's transfer energy, natom, forces, virials
    CALL mp_sum(energies, force_env%para_env%group)
    CALL mp_sum(glob_natoms, force_env%para_env%group)
    ! Transfer forces
    DO iforce_eval = 1, nforce_eval
       ALLOCATE(global_forces(iforce_eval)%forces(3,glob_natoms(iforce_eval)),stat=stat)
       CPASSERT(stat==0)
       global_forces(iforce_eval)%forces = 0.0_dp
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
       CALL mp_sum(virials(iforce_eval)%virial%pv_fock_4c, force_env%para_env%group)
       CALL mp_sum(virials(iforce_eval)%virial%pv_constraint, force_env%para_env%group)
       !Transfer results
       source = 0
       IF (ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) THEN
          IF ( force_env%sub_force_env(iforce_eval)%force_env%para_env%mepos==&
               force_env%sub_force_env(iforce_eval)%force_env%para_env%source)&
                source=force_env%para_env%mepos
       ENDIF
       CALL mp_sum(source, force_env%para_env%group)
       CALL cp_results_mp_bcast(results(iforce_eval)%results, source, force_env%para_env)
    END DO

    force_env%mixed_env%energies = energies
    ! Start combining the different sub_force_env
    CALL get_mixed_env(mixed_env=force_env%mixed_env,&
                       mixed_energy=mixed_energy)

    !NB: do this for all MIXING_TYPE values, since some need it (e.g. linear mixing 
    !NB if the first system has fewer atoms than the second)
    DO iparticle = 1, SIZE(particles_mix%els)
       particles_mix%els(iparticle)%f(:) = 0.0_dp
    END DO

    CALL section_vals_val_get(mixed_section,"MIXING_TYPE",i_val=mixing_type)
    SELECT CASE(mixing_type)
    CASE(mix_linear_combination)
       ! Support offered only 2 force_eval
       CPASSERT(nforce_eval==2)
       CALL section_vals_val_get(mixed_section,"LINEAR%LAMBDA",r_val=lambda)
       mixed_energy%pot=lambda*energies(1) + (1.0_dp-lambda)*energies(2)
       ! General Mapping of forces...
       CALL mixed_map_forces(particles_mix, virial_mix, results_mix, global_forces, virials, results,&
            lambda, 1, nforce_eval, map_index, mapping_section, .TRUE.)
       CALL mixed_map_forces(particles_mix, virial_mix, results_mix, global_forces, virials, results,&
            (1.0_dp-lambda), 2, nforce_eval, map_index, mapping_section, .FALSE.)
    CASE(mix_minimum)
       ! Support offered only 2 force_eval
       CPASSERT(nforce_eval==2)
       IF (energies(1)<energies(2)) THEN
          mixed_energy%pot=energies(1)
          CALL mixed_map_forces(particles_mix, virial_mix, results_mix, global_forces, virials, results,&
               1.0_dp, 1, nforce_eval, map_index, mapping_section, .TRUE.)
       ELSE
          mixed_energy%pot=energies(2)
          CALL mixed_map_forces(particles_mix, virial_mix, results_mix, global_forces, virials, results,&
               1.0_dp, 2, nforce_eval, map_index, mapping_section, .TRUE.)
       ENDIF
    CASE(mix_coupled)
       ! Support offered only 2 force_eval
       CPASSERT(nforce_eval==2)
       CALL section_vals_val_get(mixed_section,"COUPLING%COUPLING_PARAMETER",&
            r_val=coupling_parameter)
       sd = SQRT((energies(1)-energies(2))**2+4.0_dp*coupling_parameter**2)
       der_1=(1.0_dp-(1.0_dp/(2.0_dp*sd))*2.0_dp*(energies(1)-energies(2)))/2.0_dp
       der_2=(1.0_dp+(1.0_dp/(2.0_dp*sd))*2.0_dp*(energies(1)-energies(2)))/2.0_dp
       mixed_energy%pot=(energies(1)+energies(2)-sd)/2.0_dp
       ! General Mapping of forces...
       CALL mixed_map_forces(particles_mix, virial_mix, results_mix, global_forces, virials, results,&
            der_1, 1, nforce_eval, map_index, mapping_section, .TRUE.)
       CALL mixed_map_forces(particles_mix, virial_mix, results_mix, global_forces, virials, results,&
            der_2, 2, nforce_eval, map_index, mapping_section, .FALSE.)
    CASE(mix_restrained)
       ! Support offered only 2 force_eval
       CPASSERT(nforce_eval==2)
       CALL section_vals_val_get(mixed_section,"RESTRAINT%RESTRAINT_TARGET",&
            r_val=restraint_target)
       CALL section_vals_val_get(mixed_section,"RESTRAINT%RESTRAINT_STRENGTH",&
            r_val=restraint_strength)
       mixed_energy%pot=energies(1)+restraint_strength*(energies(1)-energies(2)-restraint_target)**2
       der_2 = -2.0_dp*restraint_strength*(energies(1)-energies(2)-restraint_target)
       der_1 = 1.0_dp - der_2
       ! General Mapping of forces...
       CALL mixed_map_forces(particles_mix, virial_mix, results_mix, global_forces, virials, results,&
            der_1, 1, nforce_eval, map_index, mapping_section, .TRUE.)
       CALL mixed_map_forces(particles_mix, virial_mix, results_mix, global_forces, virials, results,&
            der_2, 2, nforce_eval, map_index, mapping_section, .FALSE.)
    CASE(mix_generic)
       ! Support any number of force_eval sections
       gen_section => section_vals_get_subs_vals(mixed_section,"GENERIC")
       CALL get_generic_info(gen_section, "MIXING_FUNCTION", coupling_function, force_env%mixed_env%par,&
            force_env%mixed_env%val, energies)
       CALL initf(1)
       CALL parsef(1,TRIM(coupling_function),force_env%mixed_env%par)
       ! Now the hardest part.. map energy with corresponding force_eval
       mixed_energy%pot= evalf(1,force_env%mixed_env%val)
       CPASSERT(EvalErrType<=0)
       CALL zero_virial(virial_mix, reset=.FALSE.)
       CALL cp_results_erase(results_mix)
       DO iforce_eval = 1, nforce_eval
          CALL section_vals_val_get(gen_section,"DX",r_val=dx)
          CALL section_vals_val_get(gen_section,"ERROR_LIMIT",r_val=lerr)
          dedf = evalfd(1,iforce_eval,force_env%mixed_env%val,dx,err)
          IF (ABS(err)>lerr) THEN
             WRITE(this_error,"(A,G12.6,A)")"(",err,")"
             WRITE(def_error,"(A,G12.6,A)")"(",lerr,")"
             CALL compress(this_error,.TRUE.)
             CALL compress(def_error,.TRUE.)
             CALL cp_assert(.FALSE.,cp_warning_level,-300,routineP,&
                  'ASSERTION (cond) failed at line '//cp_to_string(__LINE__)//&
                  ' Error '//TRIM(this_error)//' in computing numerical derivatives larger then'//&
                  TRIM(def_error)//' .',only_ionode=.TRUE.)
          END IF
          ! General Mapping of forces...
          CALL mixed_map_forces(particles_mix, virial_mix, results_mix, global_forces, virials, results,&
               dedf, iforce_eval, nforce_eval, map_index, mapping_section, .FALSE.)
          force_env%mixed_env%val(iforce_eval) = energies(iforce_eval)
       END DO
       ! Let's store the needed information..
       force_env%mixed_env%dx  = dx
       force_env%mixed_env%lerr= lerr
       force_env%mixed_env%coupling_function = TRIM(coupling_function)
       CALL finalizef()
    CASE DEFAULT
       CPASSERT(.FALSE.)
    END SELECT
    !Simply deallocate and loose the pointer references..
    DO iforce_eval = 1, nforce_eval
       DEALLOCATE(global_forces(iforce_eval)%forces,stat=stat)
       CPASSERT(stat==0)
       CALL virial_release(virials(iforce_eval)%virial)
       CALL cp_result_release(results(iforce_eval)%results)
    END DO
    DEALLOCATE(global_forces, stat=stat)
    CPASSERT(stat==0)
    DEALLOCATE(subsystems, stat=stat)
    CPASSERT(stat==0)
    DEALLOCATE(particles, stat=stat)
    CPASSERT(stat==0)
    DEALLOCATE(energies, stat=stat)
    CPASSERT(stat==0)
    DEALLOCATE(glob_natoms, stat=stat)
    CPASSERT(stat==0)
    DEALLOCATE(virials, stat=stat)
    CPASSERT(stat==0)
    DEALLOCATE(results, stat=stat)
    CPASSERT(stat==0)
    ! Print Section
    unit_nr=cp_print_key_unit_nr(logger,mixed_section,"PRINT%DIPOLE",&
         extension=".data",middle_name="MIXED_DIPOLE",log_filename=.FALSE.)
    IF (unit_nr>0) THEN
       description ='[DIPOLE]'
       dip_exists = test_for_result(results=results_mix,description=description)
       IF (dip_exists) THEN
         CALL get_results(results=results_mix,description=description,values=dip_mix)
         WRITE(unit_nr,'(/,1X,A,T48,3F11.6)')"MIXED ENV| DIPOLE  ( A.U.)|",dip_mix
         WRITE(unit_nr,'(  1X,A,T48,3F11.6)')"MIXED ENV| DIPOLE  (Debye)|",dip_mix*debye
       ELSE
         WRITE(unit_nr,*) "NO FORCE_EVAL section calculated the dipole"
       END IF
    END IF
    CALL cp_print_key_finished_output(unit_nr,logger,mixed_section,"PRINT%DIPOLE")
  END SUBROUTINE mixed_energy_forces

END MODULE force_env_methods
