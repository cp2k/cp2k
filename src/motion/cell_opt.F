!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief performs CELL optimization
!> \par History
!>      03.2008 - Teodoro Laino [tlaino] - University of Zurich - Cell Optimization
! **************************************************************************************************
MODULE cell_opt
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                              cp_iterate,&
                                              cp_rm_iter_level
   USE force_env_types,                 ONLY: force_env_type
   USE geo_opt,                         ONLY: gopt_dispatch_optimizer
   USE global_types,                    ONLY: global_environment_type
   USE gopt_f_methods,                  ONLY: gopt_f_create_x0
   USE gopt_f_types,                    ONLY: gopt_f_create,&
                                              gopt_f_release,&
                                              gopt_f_type
   USE gopt_param_types,                ONLY: gopt_param_read,&
                                              gopt_param_type
   USE input_constants,                 ONLY: default_cell_method_id
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cell_opt'

   PUBLIC :: cp_cell_opt

CONTAINS

! **************************************************************************************************
!> \brief Main driver to perform geometry optimization
!> \param force_env ...
!> \param globenv ...
!> \author Teodoro Laino [tlaino] - University of Zurich - 03.2008
! **************************************************************************************************
   SUBROUTINE cp_cell_opt(force_env, globenv)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(global_environment_type), POINTER             :: globenv

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cp_cell_opt'

      INTEGER                                            :: handle, step_start_val
      REAL(KIND=dp), DIMENSION(:), POINTER               :: x0
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(gopt_f_type), POINTER                         :: gopt_env
      TYPE(gopt_param_type), POINTER                     :: gopt_param
      TYPE(section_vals_type), POINTER                   :: geo_section, root_section

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()
      CPASSERT(ASSOCIATED(force_env))
      CPASSERT(ASSOCIATED(globenv))
      NULLIFY (gopt_param, gopt_env, x0)
      root_section => force_env%root_section
      geo_section => section_vals_get_subs_vals(root_section, "MOTION%CELL_OPT")

      ALLOCATE (gopt_param)
      CALL gopt_param_read(gopt_param, geo_section, type_id=default_cell_method_id)
      CALL gopt_f_create(gopt_env, gopt_param, force_env=force_env, globenv=globenv, &
                         geo_opt_section=geo_section)
      CALL gopt_f_create_x0(gopt_env, x0)

      CALL section_vals_val_get(geo_section, "STEP_START_VAL", i_val=step_start_val)
      CALL cp_add_iter_level(logger%iter_info, "CELL_OPT")
      CALL cp_iterate(logger%iter_info, iter_nr=step_start_val)
      CALL gopt_dispatch_optimizer(force_env, globenv, gopt_param, gopt_env, &
                                   geo_section, x0)
      CALL cp_rm_iter_level(logger%iter_info, "CELL_OPT")

      ! Reset counter for next iteration
      CALL section_vals_val_set(geo_section, "STEP_START_VAL", i_val=0)
      DEALLOCATE (x0)
      CALL gopt_f_release(gopt_env)
      DEALLOCATE (gopt_param)
      CALL timestop(handle)

   END SUBROUTINE cp_cell_opt

END MODULE cell_opt
