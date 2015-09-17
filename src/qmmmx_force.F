!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculates QM/MM energy and forces with Force-Mixing
!> \par History
!>      2015 Factored out of force_env_methods.F
!> \author Ole Schuett
! *****************************************************************************
MODULE qmmmx_force
  USE cell_types,                      ONLY: cell_type
  USE cp_subsys_types,                 ONLY: cp_subsys_type
  USE fist_environment_types,          ONLY: fist_env_get
  USE input_constants,                 ONLY: do_fm_mom_conserv_QM,&
                                             do_fm_mom_conserv_buffer,&
                                             do_fm_mom_conserv_core,&
                                             do_fm_mom_conserv_equal_a,&
                                             do_fm_mom_conserv_equal_f,&
                                             do_fm_mom_conserv_none
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE particle_types,                  ONLY: particle_type
  USE qmmm_force,                      ONLY: qmmm_calc_energy_force
  USE qmmm_types,                      ONLY: qmmm_env_get,&
                                             qmmm_env_type
  USE qmmm_types_low,                  ONLY: force_mixing_label_QM_core,&
                                             force_mixing_label_QM_dynamics,&
                                             force_mixing_label_buffer
  USE qmmm_util,                       ONLY: apply_qmmm_unwrap,&
                                             apply_qmmm_wrap
  USE qmmmx_types,                     ONLY: qmmmx_env_type
  USE qmmmx_util,                      ONLY: apply_qmmmx_translate
  USE qs_environment_types,            ONLY: get_qs_env
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmmx_force'

  PUBLIC :: qmmmx_calc_energy_force

CONTAINS

! *****************************************************************************
!> \brief calculates the qm/mm energy and forces
!> \param qmmmx_env ...
!> \param calc_force if also the forces should be calculated
!> \param consistent_energies ...
!> \param linres ...
!> \param require_consistent_energy_force ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2004 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE qmmmx_calc_energy_force(qmmmx_env,calc_force,consistent_energies,linres,&
                                    require_consistent_energy_force,error)
    TYPE(qmmmx_env_type), POINTER            :: qmmmx_env
    LOGICAL, INTENT(IN)                      :: calc_force, &
                                                consistent_energies, linres
    LOGICAL, INTENT(IN), OPTIONAL :: require_consistent_energy_force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmmx_calc_energy_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip, mom_conserv_min_label, &
                                                mom_conserv_n, &
                                                mom_conserv_region, &
                                                mom_conserv_type
    INTEGER, POINTER                         :: cur_indices(:), cur_labels(:)
    LOGICAL                                  :: failure
    REAL(dp)                                 :: delta_a(3), delta_f(3), &
                                                mom_conserv_mass, total_f(3)
    TYPE(cp_subsys_type), POINTER            :: subsys_primary, &
                                                subsys_qmmm_core, &
                                                subsys_qmmm_extended
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_primary, &
                                                particles_qmmm_core, &
                                                particles_qmmm_extended
    TYPE(section_vals_type), POINTER         :: force_env_section

    IF (PRESENT(require_consistent_energy_force)) THEN
       CALL cp_assert(.NOT.require_consistent_energy_force,&
            cp_failure_level,cp_assertion_failed,routineP,&
            "qmmmx_energy_and_forces got require_consistent_energy_force but force mixing is active. "//&
            CPSourceFileRef,error,failure)
    ENDIF

    ! Possibly translate the system
    CALL apply_qmmmx_translate(qmmmx_env, error)

    ! actual energy force calculation
    CALL qmmmx_calc_energy_force_low(qmmmx_env%ext, calc_force, consistent_energies, linres, "ext", error=error)
    CALL qmmmx_calc_energy_force_low(qmmmx_env%core, calc_force, consistent_energies, linres, "core", error=error)

    ! get forces from subsys of each sub force env
    CALL qmmm_env_get(qmmmx_env%core, subsys=subsys_qmmm_core, error=error)
    CALL qmmm_env_get(qmmmx_env%ext, subsys=subsys_qmmm_extended, error=error)

    CALL get_qs_env(qmmmx_env%ext%qs_env, input=force_env_section, error=error)
    CALL section_vals_val_get(force_env_section,"QMMM%FORCE_MIXING%RESTART_INFO%INDICES",i_vals=cur_indices,error=error)
    CALL section_vals_val_get(force_env_section,"QMMM%FORCE_MIXING%RESTART_INFO%LABELS",i_vals=cur_labels,error=error)

    particles_qmmm_extended => subsys_qmmm_extended%particles%els
    particles_qmmm_core => subsys_qmmm_core%particles%els
    DO ip=1,SIZE(cur_indices)
       IF (cur_labels(ip) >= force_mixing_label_QM_dynamics) THEN ! this is a QM atom
         ! copy (QM) force from extended calculation
         particles_qmmm_core(cur_indices(ip))%f=particles_qmmm_extended(cur_indices(ip))%f
       ENDIF
    END DO

    ! zero momentum
    CALL section_vals_val_get(force_env_section,"QMMM%FORCE_MIXING%MOMENTUM_CONSERVATION_TYPE",&
                              i_val=mom_conserv_type,error=error)
    IF (mom_conserv_type /= do_fm_mom_conserv_none) THEN
       CALL section_vals_val_get(force_env_section,"QMMM%FORCE_MIXING%MOMENTUM_CONSERVATION_REGION",&
                                 i_val=mom_conserv_region,error=error)

       IF (mom_conserv_region == do_fm_mom_conserv_core) THEN
          mom_conserv_min_label = force_mixing_label_QM_core
       ELSEIF (mom_conserv_region == do_fm_mom_conserv_QM) THEN
          mom_conserv_min_label = force_mixing_label_QM_dynamics
       ELSEIF (mom_conserv_region == do_fm_mom_conserv_buffer) THEN
          mom_conserv_min_label = force_mixing_label_buffer
       ELSE
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"Got unknown MOMENTUM_CONSERVATION_REGION (not CORE, QM, or BUFFER) !"//&
             CPSourceFileRef,&
            error,failure)
       ENDIF

       total_f = 0.0_dp
       DO ip=1, SIZE(particles_qmmm_core)
         total_f(1:3) = total_f(1:3) + particles_qmmm_core(ip)%f(1:3)
       END DO
       IF (mom_conserv_type == do_fm_mom_conserv_equal_f) THEN
          mom_conserv_n = COUNT(cur_labels >= mom_conserv_min_label)
          delta_f = total_f/mom_conserv_n
          DO ip=1, SIZE(cur_indices)
            IF (cur_labels(ip) >= mom_conserv_min_label) THEN
              particles_qmmm_core(cur_indices(ip))%f = particles_qmmm_core(cur_indices(ip))%f - delta_f
            ENDIF
          END DO
       ELSE IF (mom_conserv_type == do_fm_mom_conserv_equal_a) THEN
          mom_conserv_mass = 0.0_dp
          DO ip=1, SIZE(cur_indices)
           IF (cur_labels(ip) >= mom_conserv_min_label) &
             mom_conserv_mass = mom_conserv_mass + particles_qmmm_core(cur_indices(ip))%atomic_kind%mass
          END DO
          delta_a = total_f/mom_conserv_mass
          DO ip=1, SIZE(cur_indices)
            IF (cur_labels(ip) >= mom_conserv_min_label) THEN
              particles_qmmm_core(cur_indices(ip))%f = particles_qmmm_core(cur_indices(ip))%f - &
                 particles_qmmm_core(cur_indices(ip))%atomic_kind%mass * delta_a
            ENDIF
          END DO
       ENDIF
    ENDIF

    CALL qmmm_env_get(qmmmx_env%ext, subsys=subsys_primary, error=error)
    particles_primary => subsys_primary%particles%els
    DO ip=1,SIZE(particles_qmmm_core)
       particles_primary(ip)%f=particles_qmmm_core(ip)%f
    END DO

  END SUBROUTINE qmmmx_calc_energy_force

! *****************************************************************************
!> \brief ...
!> \param qmmm_env ...
!> \param calc_force ...
!> \param consistent_energies ...
!> \param linres ...
!> \param label ...
!> \param error ...
! *****************************************************************************
 SUBROUTINE qmmmx_calc_energy_force_low(qmmm_env,calc_force,consistent_energies,linres,label,error)
    TYPE(qmmm_env_type), POINTER             :: qmmm_env
    LOGICAL, INTENT(IN)                      :: calc_force, &
                                                consistent_energies, linres
    CHARACTER(*)                             :: label
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmmx_calc_energy_force_low', &
      routineP = moduleN//':'//routineN

    CHARACTER(default_string_length)         :: new_restart_fn, &
                                                new_restart_hist_fn, &
                                                old_restart_fn, &
                                                old_restart_hist_fn
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index
    LOGICAL                                  :: saved_do_translate
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: saved_pos
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_subsys_type), POINTER            :: subsys_mm, subsys_qm
    TYPE(section_vals_type), POINTER         :: force_env_section

    NULLIFY(mm_cell,subsys_qm,subsys_mm,qm_atom_index)

    CALL get_qs_env(qmmm_env%qs_env, input=force_env_section, error=error)

    ! rewrite RESTART%FILENAME
    CALL section_vals_val_get(force_env_section, "DFT%SCF%PRINT%RESTART%FILENAME",&
                              c_val=old_restart_fn, error=error)
    new_restart_fn = TRIM(old_restart_fn)//"-"//TRIM(label)
    CALL section_vals_val_set(force_env_section, "DFT%SCF%PRINT%RESTART%FILENAME",&
                              c_val=new_restart_fn, error=error)

    ! rewrite RESTART_HISTORY%FILENAME
    CALL section_vals_val_get(force_env_section, "DFT%SCF%PRINT%RESTART_HISTORY%FILENAME",&
                              c_val=old_restart_hist_fn, error=error)
    new_restart_hist_fn = TRIM(old_restart_hist_fn)//"-"//TRIM(label)
    CALL section_vals_val_set(force_env_section, "DFT%SCF%PRINT%RESTART_HISTORY%FILENAME",&
                              c_val=new_restart_hist_fn, error=error)

    ! wrap positions before QM/MM calculation.
    ! Required if diffusion causes atoms outside of periodic box get added to QM
    CALL fist_env_get(qmmm_env%fist_env, cell=mm_cell,subsys=subsys_mm,error=error)
    CALL get_qs_env(qmmm_env%qs_env, cp_subsys=subsys_qm,error=error)
    qm_atom_index => qmmm_env%qm%qm_atom_index
    CALL apply_qmmm_wrap(subsys_mm, mm_cell, subsys_qm, qm_atom_index, saved_pos, error)

    ! Turn off box translation, it was already performed by apply_qmmmx_translate(),
    ! the particles coordinates will still be copied from MM to QM.
    saved_do_translate = qmmm_env%qm%do_translate
    qmmm_env%qm%do_translate = .FALSE.

    ! actual energy force calculation
    CALL qmmm_calc_energy_force(qmmm_env, calc_force, consistent_energies, linres, error)

    ! restore do_translate
    qmmm_env%qm%do_translate = saved_do_translate

    ! restore unwrapped positions
    CALL apply_qmmm_unwrap(subsys_mm, subsys_qm, qm_atom_index, saved_pos, error)

    ! restore RESTART filenames
    CALL section_vals_val_set(force_env_section, "DFT%SCF%PRINT%RESTART%FILENAME",&
                              c_val=old_restart_fn, error=error)
    CALL section_vals_val_set(force_env_section, "DFT%SCF%PRINT%RESTART_HISTORY%FILENAME",&
                              c_val=old_restart_hist_fn, error=error)

  END SUBROUTINE qmmmx_calc_energy_force_low

END MODULE qmmmx_force
