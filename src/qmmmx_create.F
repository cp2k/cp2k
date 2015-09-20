!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Initialize a QM/MM calculation with Force-Mixing
!> \author Ole Schuett
! *****************************************************************************
MODULE qmmmx_create
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_release,&
                                             section_vals_type
  USE qmmm_create,                     ONLY: qmmm_env_create
  USE qmmm_types,                      ONLY: qmmm_env_get,&
                                             qmmm_env_release,&
                                             qmmm_env_type
  USE qmmmx_types,                     ONLY: qmmmx_env_type
  USE qmmmx_util,                      ONLY: setup_force_mixing_qmmm_sections,&
                                             update_force_mixing_labels
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmmx_create'

  PUBLIC :: qmmmx_env_create

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param qmmmx_env ...
!> \param root_section ...
!> \param para_env ...
!> \param globenv ...
!> \param force_env_section ...
!> \param subsys_section ...
!> \param use_motion_section ...
!> \par History
!>      02.2012 created [noam]
!> \author Noam Bernstein
! *****************************************************************************
  SUBROUTINE  qmmmx_env_create(qmmmx_env, root_section, para_env, globenv,&
       force_env_section, subsys_section, use_motion_section)
    TYPE(qmmmx_env_type), POINTER            :: qmmmx_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                subsys_section
    LOGICAL, INTENT(IN)                      :: use_motion_section

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmmx_env_create', &
      routineP = moduleN//':'//routineN

    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(qmmm_env_type), POINTER             :: dummy_qmmm_env
    TYPE(section_vals_type), POINTER         :: qmmm_core_section, &
                                                qmmm_extended_section, &
                                                qmmm_section

    NULLIFY(dummy_qmmm_env)

    qmmm_section => section_vals_get_subs_vals(force_env_section,"QMMM")

    CALL qmmm_env_create(dummy_qmmm_env, root_section, para_env, globenv,&
                         force_env_section, qmmm_section, subsys_section, use_motion_section, &
                         ignore_outside_box = .TRUE.)
    CALL qmmm_env_get(dummy_qmmm_env, subsys=subsys)

    CALL update_force_mixing_labels(subsys, qmmm_section)

    ! using CUR_INDICES and CUR_LABELS, create appropriate QM_KIND sections for two QM/MM calculations
    CALL setup_force_mixing_qmmm_sections(subsys, qmmm_section, qmmm_core_section, qmmm_extended_section)


    ALLOCATE(qmmmx_env)
    CALL qmmm_env_create(qmmmx_env%core, root_section, para_env, globenv,&
         force_env_section, qmmm_core_section, subsys_section, use_motion_section, &
         ignore_outside_box = .TRUE.)

    CALL qmmm_env_create(qmmmx_env%ext, root_section, para_env, globenv,&
         force_env_section, qmmm_extended_section, subsys_section, use_motion_section, &
         ignore_outside_box = .TRUE.)

    CALL section_vals_release(qmmm_core_section)
    CALL section_vals_release(qmmm_extended_section)
    CALL qmmm_env_release(dummy_qmmm_env)

  END SUBROUTINE qmmmx_env_create

END MODULE qmmmx_create
