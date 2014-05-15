!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Only ep_f77 is supposed to use these methods
!> \author Ole Schuett
! *****************************************************************************

! *****************************************************************************
  SUBROUTINE ep_f77_create_low()
    USE ep_types,                        ONLY: ep_env_type
    USE ep_methods,                      ONLY: ep_env_create
    USE ep_f77,                          ONLY: param_ep_env_id,&
                                               param_root_section,&
                                               param_input_declaration,&
                                               param_para_env,&
                                               param_globenv,&
                                               param_error
   IMPLICIT NONE
    TYPE(ep_env_type), POINTER               :: ep_env
    NULLIFY(ep_env)
    CALL ep_env_create(ep_env,param_root_section,param_input_declaration,&
           param_para_env, param_globenv, param_error)
    param_ep_env_id = ep_env%id_nr
  END SUBROUTINE ep_f77_create_low


! *****************************************************************************
  SUBROUTINE ep_f77_get_energy_low()
    USE ep_f77,                          ONLY: param_ep_env_id,&
                                               param_energy
    USE ep_types,                        ONLY: ep_env_type
    USE ep_methods,                      ONLY: ep_envs_get_ep_env
    IMPLICIT NONE
    TYPE(ep_env_type), POINTER               :: ep_env
    ep_env => ep_envs_get_ep_env(param_ep_env_id)
    ALLOCATE(param_energy)
    param_energy = ep_env%energy%e_tot
  END SUBROUTINE ep_f77_get_energy_low


! *****************************************************************************
  SUBROUTINE ep_f77_get_main_qs_env_low()
    USE ep_f77,                          ONLY: param_ep_env_id,&
                                               param_main_qs_env
    USE ep_types,                        ONLY: ep_env_type
    USE ep_methods,                      ONLY: ep_envs_get_ep_env
    IMPLICIT NONE
    TYPE(ep_env_type), POINTER               :: ep_env
    ep_env => ep_envs_get_ep_env(param_ep_env_id)
    param_main_qs_env => ep_env%main_qs_env
  END SUBROUTINE ep_f77_get_main_qs_env_low


! *****************************************************************************
  SUBROUTINE ep_f77_retain_low()
    USE ep_f77,                          ONLY: param_ep_env_id,&
                                               param_error
    USE ep_types,                        ONLY: ep_env_type
    USE ep_methods,                      ONLY: ep_envs_get_ep_env,&
                                               ep_env_retain
    IMPLICIT NONE
    TYPE(ep_env_type), POINTER               :: ep_env
    ep_env => ep_envs_get_ep_env(param_ep_env_id)
    CALL ep_env_retain(ep_env, param_error)
  END SUBROUTINE ep_f77_retain_low


! *****************************************************************************
  SUBROUTINE ep_f77_release_low()
    USE ep_f77,                          ONLY: param_ep_env_id,&
                                               param_error
    USE ep_types,                        ONLY: ep_env_type
    USE ep_methods,                      ONLY: ep_envs_get_ep_env,&
                                               ep_env_release
    IMPLICIT NONE
    TYPE(ep_env_type), POINTER               :: ep_env
    INTEGER                                  :: ref_count
    ep_env => ep_envs_get_ep_env(param_ep_env_id)
    ref_count = ep_env%ref_count
    CALL ep_env_release(ep_env, param_error)
    IF(ref_count==1)&
       DEALLOCATE(param_ep_env_id)
  END SUBROUTINE ep_f77_release_low


! *****************************************************************************
  SUBROUTINE ep_f77_calc_e_f_low()
    USE ep_f77,                          ONLY: param_ep_env_id,&
                                               param_calc_f,&
                                               param_error
    USE ep_types,                        ONLY: ep_env_type
    USE ep_methods,                      ONLY: ep_envs_get_ep_env,&
                                               ep_env_calc_e_f
    IMPLICIT NONE
    TYPE(ep_env_type), POINTER               :: ep_env
    ep_env => ep_envs_get_ep_env(param_ep_env_id)
    CALL ep_env_calc_e_f(ep_env,param_calc_f,param_error)
  END SUBROUTINE ep_f77_calc_e_f_low


!EOF
