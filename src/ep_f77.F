!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Ugly f77-hack to resolve the following circular dependency:
!>           force_env -> ep_method -> replica_env -> force_env
!> \author Ole Schuett
! *****************************************************************************
MODULE ep_f77
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_type,&
                                             section_vals_type
  USE iso_c_binding
  USE kinds,                           ONLY: dp
  USE qs_environment_types,            ONLY: qs_environment_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ep_f77_calc_e_f, ep_f77_create, ep_f77_release, ep_f77_retain
  PUBLIC :: ep_f77_get_main_qs_env, ep_f77_get_energy

  ! public global variables used to hand over arguments
  ! only ep_f77_low is supposed to use them
  INTEGER, POINTER, SAVE                   :: param_ep_env_id => Null()
  TYPE(section_vals_type), POINTER, SAVE   :: param_root_section => Null()
  TYPE(section_type), POINTER, SAVE        :: param_input_declaration => Null()
  TYPE(cp_para_env_type), POINTER, SAVE    :: param_para_env => Null()
  TYPE(global_environment_type), POINTER, SAVE :: param_globenv => Null()
  TYPE(cp_error_type), POINTER, SAVE       :: param_error => Null()
  LOGICAL, SAVE                            :: param_calc_f = .FALSE.
  TYPE(qs_environment_type), POINTER, SAVE :: param_main_qs_env => Null()
  REAL(KIND=dp), POINTER, SAVE             :: param_energy

  PUBLIC:: param_ep_env_id, param_root_section, param_input_declaration
  PUBLIC:: param_para_env, param_globenv, param_error, param_calc_f
  PUBLIC:: param_main_qs_env, param_energy

CONTAINS

! *****************************************************************************
!> \brief wrapper for ep_env_create()
!> \param ep_env_id ...
!> \param root_section ...
!> \param input_declaration ...
!> \param para_env ...
!> \param globenv ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE ep_f77_create(ep_env_id,root_section,input_declaration,para_env,globenv, error)
    INTEGER, POINTER                         :: ep_env_id
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(section_type), POINTER              :: input_declaration
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      TARGET                                 :: error

    param_root_section       => root_section
    param_input_declaration  => input_declaration
    param_para_env           => para_env
    param_globenv            => globenv
    param_error              => error
    ALLOCATE(param_ep_env_id)

    CALL ep_f77_create_low()

    ALLOCATE(ep_env_id)
    ep_env_id = param_ep_env_id

    DEALLOCATE(param_ep_env_id)
    NULLIFY(param_root_section, param_input_declaration)
    NULLIFY(param_para_env, param_globenv, param_error)
  END SUBROUTINE ep_f77_create


! *****************************************************************************
!> \brief returns ep_env%energy%e_tot
!> \param ep_env_id ...
!> \retval energy ...
!> \author Ole Schuett
! *****************************************************************************
  FUNCTION ep_f77_get_energy(ep_env_id) RESULT(energy)
    INTEGER, POINTER                         :: ep_env_id
    REAL(KIND=dp)                            :: energy

    param_ep_env_id => ep_env_id
    CALL ep_f77_get_energy_low()
    energy = param_energy
    DEALLOCATE(param_energy)
    NULLIFY(param_ep_env_id)
  END FUNCTION ep_f77_get_energy


! *****************************************************************************
!> \brief returns ep_env%main_qs_env
!> \param ep_env_id ...
!> \param main_qs_env ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE ep_f77_get_main_qs_env(ep_env_id, main_qs_env)
    INTEGER, POINTER                         :: ep_env_id
    TYPE(qs_environment_type), POINTER       :: main_qs_env

    param_ep_env_id => ep_env_id

    CALL ep_f77_get_main_qs_env_low()

    main_qs_env => param_main_qs_env
    NULLIFY(param_ep_env_id, param_main_qs_env)
  END SUBROUTINE ep_f77_get_main_qs_env


! *****************************************************************************
!> \brief wrapper for ep_env_retain()
!> \param ep_env_id ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE ep_f77_retain(ep_env_id, error)
    INTEGER, POINTER                         :: ep_env_id
    TYPE(cp_error_type), INTENT(inout), &
      TARGET                                 :: error

    param_ep_env_id => ep_env_id
    param_error     => error
    CALL ep_f77_retain_low()
    NULLIFY(param_ep_env_id, param_error)
  END SUBROUTINE ep_f77_retain


! *****************************************************************************
!> \brief wrapper for ep_env_release()
!> \param ep_env_id ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE ep_f77_release(ep_env_id, error)
    INTEGER, POINTER                         :: ep_env_id
    TYPE(cp_error_type), INTENT(inout), &
      TARGET                                 :: error

    param_ep_env_id => ep_env_id
    param_error     => error
    CALL ep_f77_release_low()
    NULLIFY(param_ep_env_id, param_error)
  END SUBROUTINE ep_f77_release


! *****************************************************************************
!> \brief wrapper for ep_env_calc_e_f()
!> \param ep_env_id ...
!> \param calc_f ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE ep_f77_calc_e_f(ep_env_id, calc_f, error)
    INTEGER, POINTER                         :: ep_env_id
    LOGICAL, INTENT(in)                      :: calc_f
    TYPE(cp_error_type), INTENT(inout), &
      TARGET                                 :: error

    param_ep_env_id          => ep_env_id
    param_calc_f             = calc_f
    param_error              => error

    CALL ep_f77_calc_e_f_low()
    NULLIFY(param_ep_env_id, param_error)
  END SUBROUTINE ep_f77_calc_e_f

END MODULE ep_f77
