!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Initialize a qs_env for kpoint calculations starting from a gamma point qs_env
!> \par History
!>      11.2016 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE qs_gamma2kp
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_subsys_types,                 ONLY: cp_subsys_type
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_create,&
                                              kpoint_type
   USE qs_energy_init,                  ONLY: qs_energies_init
   USE qs_environment,                  ONLY: qs_init
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_env_create,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
   USE qs_scf_initialization,           ONLY: qs_scf_env_init_basic
   USE qs_scf_types,                    ONLY: qs_scf_env_type,&
                                              scf_env_release
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_gamma2kp'

   PUBLIC :: create_kp_from_gamma

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param qs_env_kp ...
! **************************************************************************************************
   SUBROUTINE create_kp_from_gamma(qs_env, qs_env_kp)
      TYPE(qs_environment_type), POINTER                 :: qs_env, qs_env_kp

      CHARACTER(len=*), PARAMETER :: routineN = 'create_kp_from_gamma', &
         routineP = moduleN//':'//routineN

      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_subsys_type), POINTER                      :: cp_subsys
      TYPE(kpoint_type), POINTER                         :: kpoint
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(section_vals_type), POINTER                   :: force_env_section, subsys_section

      !--------------------------------------------------------------------------------------------!

      CPASSERT(.NOT. ASSOCIATED(qs_env_kp))

      CALL get_qs_env(qs_env, &
                      para_env=para_env, &
                      input=force_env_section, &
                      cp_subsys=cp_subsys)

      NULLIFY (subsys_section)

      NULLIFY (kpoint)
      CALL kpoint_create(kpoint)
      kpoint%kp_scheme = "GAMMA"
      kpoint%symmetry = .FALSE.
      kpoint%verbose = .FALSE.
      kpoint%full_grid = .TRUE.
      kpoint%eps_geo = 1.0e-6_dp
      kpoint%use_real_wfn = .TRUE.
      kpoint%parallel_group_size = 0

      CALL qs_env_create(qs_env_kp)
      CALL qs_init(qs_env_kp, para_env, cp_subsys=cp_subsys, kpoint_env=kpoint, &
                   force_env_section=force_env_section, subsys_section=subsys_section, &
                   use_motion_section=.FALSE.)

      CALL qs_energies_init(qs_env_kp, calc_forces=.FALSE.)

      NULLIFY (scf_env)
      CALL qs_scf_env_init_basic(qs_env_kp, scf_env)
      CALL set_qs_env(qs_env_kp, scf_env=scf_env)
      CALL scf_env_release(scf_env)

      CALL qs_ks_update_qs_env(qs_env_kp, print_active=.FALSE.)

   END SUBROUTINE create_kp_from_gamma

END MODULE qs_gamma2kp
